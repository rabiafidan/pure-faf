# Rabia updated this script on 21/10/2025 to make it accept statistical analysis pval threshold parameter as argument.
# threads argument added on 22/10/2025
# for pval, use the scientific notation without ^, e.g., 1e-6
# usage: Rscript MicroSEC_parameterised.R [working/output directory] [sample information tsv file] [progress bar Y/N] [pval threshold] [threads]
# MicroSEC: Microhomology-induced chimeric read-originating sequence error cleaning pipeline for FFPE samples
#
# Author: "Masachika Ikegami"
#
# ---------------------------
# CHANGELOG (25/10/2025)
# ---------------------------
# - Final .tsv.gz is written directly into wd; only wd/tmp/ is created for intermediates.
# - Streamed samtools view|sort (no on-disk SAM); CRAM uses -T <reference.fa>.
# - Bound samtools sort RAM with -m and temp location with -T + TMPDIR=wd/tmp.
# - Always (re)index slim BAM; strict failure on non-zero system() exit.
# - Fix: Filter df_mutation to BAM contigs to avoid fun_read_check() '[[1]] out of bounds'.
# - Normalize columns; fallback for adapters; safe zero-variant short-circuit.

suppressPackageStartupMessages({
  library(MicroSEC)
  library(dplyr)
  library(readr)
  library(stringr)
  library(Rsamtools)
  library(BiocGenerics)
  library(Biostrings)
  library(tools)
})

# -------- args --------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript MicroSEC_parameterised.R [wd] [sample_info.tsv] [progress Y/N] [threshold_p] [threads]")
}
wd           <- args[1]
sample_list  <- args[2]
progress_bar <- args[3]
threshold_p  <- as.numeric(args[4])
threads      <- as.numeric(args[5])
memory_per_thread <- paste0(as.character(args[6]),"M")  # e.g., "2000M" for 2GB

wd <- normalizePath(wd, mustWork = FALSE)
if (!dir.exists(wd)) stop(sprintf("Working/output directory does not exist: %s", wd))
setwd(wd)

# Read sample info
sample_info <- read.csv(sample_list, header = FALSE, stringsAsFactors = FALSE, sep = "\t")

# Init containers
msec <- NULL
homology_search <- NULL
mut_depth <- list(NULL, NULL, NULL)

# Helpers
run_cmd <- function(cmd) {
  status <- system(cmd)
  if (status != 0) stop(sprintf("Command failed (exit %s): %s", status, cmd))
}
q <- function(x) shQuote(normalizePath(x, mustWork = FALSE))

# ----------------------------- per-sample loop -----------------------------
for (sample in seq_len(nrow(sample_info))) {
  sample_name   <- sample_info[sample, 1]
  mutation_file <- sample_info[sample, 2]
  bam_file_in   <- sample_info[sample, 3]
  read_length   <- as.integer(sample_info[sample, 4])
  adapter_1     <- sample_info[sample, 5]

  # Optional columns handling
  if (sample_info[sample, 6] %in% c("Human","Mouse","hg19","hg38","mm10")) {
    adapter_2 <- adapter_1
    organism  <- sample_info[sample, 6]
    panel     <- sample_info[sample, 7]
    if (ncol(sample_info) == 8)  reference_genome <- sample_info[sample, 8]
    if (ncol(sample_info) == 9){ reference_genome <- sample_info[sample, 8]; simple_repeat_list <- sample_info[sample, 9] }
  } else {
    adapter_2 <- sample_info[sample, 6]
    organism  <- sample_info[sample, 7]
    panel     <- sample_info[sample, 8]
    if (ncol(sample_info) == 9)  reference_genome <- sample_info[sample, 9]
    if (ncol(sample_info) == 10){ reference_genome <- sample_info[sample, 9]; simple_repeat_list <- sample_info[sample,10] }
  }
  if (is.na(adapter_2) || adapter_2 == "") adapter_2 <- adapter_1

  # Paths: outputs go directly to wd/, only tmp/ is created
  tmp_dir <- file.path(wd, "tmp")
  if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  bam_file_slim     <- file.path(tmp_dir, paste0(sample_name, ".SLIM.bam"))
  bam_file_slim_bai <- paste0(bam_file_slim, ".bai")
  bed_path          <- file.path(tmp_dir, paste0(sample_name, ".regions.bed"))

  # Genome & chromosome naming
  ref_genome <- fun_load_genome(organism)
  chr_no     <- fun_load_chr_no(organism)
  if (ref_genome@user_seqnames[[1]] == "chr1") {
    chromosomes <- paste0("chr", c(seq_len(chr_no - 2), "X", "Y"))
  } else {
    chromosomes <- paste0("", c(seq_len(chr_no - 2), "X", "Y"))
  }

  # Load mutations
  df_mutation <- fun_load_mutation(mutation_file, sample_name, ref_genome, chr_no)

  # Normalize required columns & fallbacks
  req_cols <- c("Sample","Mut_type","Chr","Pos","Ref","Alt","SimpleRepeat_TRF","Neighborhood_sequence")
  missing <- setdiff(req_cols, colnames(df_mutation))
  if (length(missing) > 0) {
    # create missing columns with defaults
    for (mc in missing) df_mutation[[mc]] <- "-"
  }
  df_mutation <- df_mutation %>%
    mutate(
      Sample = as.character(Sample),
      Mut_type = as.character(Mut_type),
      Chr = as.character(Chr),
      Pos = as.integer(Pos),
      Ref = as.character(Ref),
      Alt = as.character(Alt),
      SimpleRepeat_TRF = ifelse(is.na(SimpleRepeat_TRF) | SimpleRepeat_TRF == "", "-", as.character(SimpleRepeat_TRF)),
      Neighborhood_sequence = ifelse(is.na(Neighborhood_sequence) | Neighborhood_sequence == "", "-", as.character(Neighborhood_sequence))
    ) %>%
    filter(!is.na(Pos), Pos > 0)

  # Build BED around mutations (±200; merge if next within 400bp)
  df_mutation_save <- df_mutation
  download_region <- data.frame(matrix(rep(NA,3), nrow=1))[numeric(0),]
  colnames(download_region) <- c("chrom","chromStart","chromEnd")
  for (i in chromosomes) {
    df_chr <- df_mutation_save[df_mutation_save$Chr == i,]
    if (nrow(df_chr) == 0) next
    start <- max(1, df_chr$Pos[1] - 200)
    for (k in seq_len(nrow(df_chr))) {
      if (k == nrow(df_chr) || df_chr$Pos[k+1] - df_chr$Pos[k] > 400) {
        download_region <- rbind(download_region,
                                 c(i, start - 1, df_chr$Pos[k] + 200))
        if (k < nrow(df_chr)) start <- max(1, df_chr$Pos[k+1] - 200)
      }
    }
  }
  write_tsv(download_region, bed_path, col_names = FALSE, progress = FALSE)

  # ---- samtools pipeline (streamed) ----
  mem_per_thread <- memory_per_thread
  sort_prefix <- file.path(tmp_dir, "samtools_sort")
  input_ext <- tolower(file_ext(bam_file_in))

  if (input_ext == "bam") {
    cmd <- sprintf(
      paste(
        "bash -lc 'set -euo pipefail;",
        "export TMPDIR=%s;",
        "samtools view -@ %d -b -h --no-PG %s -L %s |",
        "samtools sort -@ %d -m %s -O BAM -T %s -o %s -'",
        sep=" "
      ),
      q(tmp_dir), threads, q(bam_file_in), q(bed_path),
      threads, mem_per_thread, q(sort_prefix), q(bam_file_slim)
    )
  } else if (input_ext == "cram") {
    if (!exists("reference_genome"))
      stop("CRAM input requires reference_genome in sample info TSV.")
    cmd <- sprintf(
      paste(
        "bash -lc 'set -euo pipefail;",
        "export TMPDIR=%s;",
        "samtools view -@ %d -b -h --no-PG %s -T %s -L %s |",
        "samtools sort -@ %d -m %s -O BAM -T %s -o %s -'",
        sep=" "
      ),
      q(tmp_dir), threads, q(bam_file_in), q(reference_genome), q(bed_path),
      threads, mem_per_thread, q(sort_prefix), q(bam_file_slim)
    )
  } else {
    stop(sprintf("Unsupported input extension: %s", input_ext))
  }

  message("Trimming and sorting reads into slim BAM (tmp/)...")
  run_cmd(cmd)
  run_cmd(sprintf("bash -lc 'set -euo pipefail; samtools index %s'", q(bam_file_slim)))

  # ---- load slim BAM ----
  df_bam <- fun_load_bam(bam_file_slim)

  # ---- align mutation contigs to BAM targets (fixes [[1]] out-of-bounds) ----
  targets <- names(Rsamtools::scanBamHeader(bam_file_slim)[[1]]$targets)
  missing_chr <- setdiff(unique(df_mutation$Chr), targets)
  if (length(missing_chr) > 0)
    message(sprintf("[%s] Skipping contigs not in BAM: %s",
                  sample_name, paste(missing_chr, collapse = ", ")))
df_mutation <- dplyr::filter(df_mutation, Chr %in% targets)

  # If nothing left for this sample after subsetting, emit empty result and continue
  if (nrow(df_mutation) == 0) {
    message(sprintf("[%s] No mutations left on BAM contigs after subsetting; writing empty result.", sample_name))
    empty_path <- file.path(wd, paste0(sample_name, ".tsv.gz"))
    # Write an empty TSV with header matching fun_save’s output schema
    # Minimal safe fallback: let MicroSEC save an empty frame
    empty_msec <- data.frame()
    con <- gzfile(empty_path, "w")
    writeLines("", con) # create empty gz
    close(con)
    next
  }

  # ---- MicroSEC core ----
  result <- fun_read_check(
    df_mutation = df_mutation,
    df_bam = df_bam,
    ref_genome = ref_genome,
    sample_name = sample_name,
    read_length = read_length,
    adapter_1 = adapter_1,
    adapter_2 = adapter_2,
    short_homology_search_length = 4,
    min_homology_search = 40,
    progress_bar = progress_bar
  )

  msec            <- rbind(msec, result[[1]])
  homology_search <- rbind(homology_search, result[[2]])
  mut_depth <- list(
    rbind(mut_depth[[1]], result[[3]][[1]]),
    rbind(mut_depth[[2]], result[[3]][[2]]),
    rbind(mut_depth[[3]], result[[3]][[3]])
  )
} # end for samples

# -----------------------
# Post-processing & stats
# -----------------------
msec <- fun_homology(
  msec,
  homology_search,
  min_homology_search = 40,
  ref_genome,
  chr_no,
  progress_bar = progress_bar
)
msec <- fun_summary(msec)
msec <- fun_analysis(
  msec,
  mut_depth,
  short_homology_search_length = 4,
  min_homology_search = 40,
  threshold_p = threshold_p,
  threshold_hairpin_ratio = 0.50,
  threshold_short_length = 0.75,
  threshold_distant_homology = 0.15,
  threshold_soft_clip_ratio = 0.50,
  threshold_low_quality_rate = 0.1,
  homopolymer_length = 15
)

# ---- save final result directly to wd ----
final_out_path <- file.path(wd, paste0(sample_info[1,1], ".tsv.gz"))
fun_save(msec, final_out_path)

message("MicroSEC finished successfully.")
