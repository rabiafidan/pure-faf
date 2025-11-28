# MicroSEC_snakemake.R
#
# MicroSEC: Microhomology-induced chimeric read-originating sequence
#           error cleaning pipeline for FFPE samples
#
# Snakemake script version:
#   - Uses snakemake@input, snakemake@output, snakemake@params, snakemake@threads.
#   - Assumes a single sample in sample_info TSV.
#   - Writes final .tsv.gz exactly to snakemake@output[["msec"]].
#   - Creates tmp/ next to that output for intermediates.
#   - No setwd(), everything is absolute paths.
#
# Expected Snakemake interface:
#   input:
#     sample_info:  TSV with 1 row:
#       V1: sample_name
#       V2: mutation_file
#       V3: bam_or_cram
#       V4: read_length
#       V5: adapter_1
#       V6+: optional organism/adapter_2/panel/reference_genome/simple_repeat_list (same logic as original)
#   output:
#     msec: final MicroSEC result .tsv.gz
#   params:
#     threshold_p:     numeric p-value threshold
#     mem_per_thread: memory per thread in MB (e.g., 2000)
#   threads:
#     integer, passed to samtools and MicroSEC

log <- file(snakemake@log[[1]], open = "at")
sink(log, append = TRUE)
sink(log, type = "message", append = TRUE)
message("-----------------------------------------------")
message("-----------------------------------------------")
message(date())
##############################################
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

# ----------------- Snakemake inputs/outputs/params -----------------

# final output path (used as-is, absolute)
out_path <- normalizePath(snakemake@output[["msec"]], mustWork = FALSE)

# sample_info TSV
sample_list <- normalizePath(snakemake@input[["sample_info"]], mustWork = TRUE)

# params
progress_bar <- "Y"
threshold_p  <- as.numeric(snakemake@params[["threshold_p"]])
threads      <- as.integer(snakemake@threads)

# memory per thread (MB) -> samtools -m string
params_list <- snakemake@params
if ("mem_per_thread" %in% names(params_list)) {
  memory_per_thread <- paste0(as.character(params_list[["mem_per_thread"]]), "M")
} else {
  memory_per_thread <- "4000M"  # default 4GB if not set
}

# "wd" is only used for tmp placement next to the final output
wd <- dirname(out_path)
if (!dir.exists(wd)) {
  stop(sprintf("Working/output directory does not exist: %s", wd))
}

# ----------------- Helpers -----------------

run_cmd <- function(cmd) {
  status <- system(cmd)
  if (status != 0) stop(sprintf("Command failed (exit %s): %s", status, cmd))
}

qpath <- function(x) {
  # Normalise but do not require existence (for new files)
  shQuote(normalizePath(x, mustWork = FALSE))
}

# ----------------- Read and unpack sample_info (single sample) -----------------

sample_info <- read.csv(sample_list, header = FALSE, stringsAsFactors = FALSE, sep = "\t")
if (nrow(sample_info) != 1) {
  stop(sprintf("This script expects exactly one sample in sample_info, but found %d.", nrow(sample_info)))
}

sample_name   <- sample_info[1, 1]
mutation_file <- sample_info[1, 2]
bam_file_in   <- sample_info[1, 3]
read_length   <- as.integer(sample_info[1, 4])
adapter_1     <- sample_info[1, 5]

# Optional columns handling (same logic as original)
if (sample_info[1, 6] %in% c("Human","Mouse","hg19","hg38","mm10")) {
  adapter_2 <- adapter_1
  organism  <- sample_info[1, 6]
  panel     <- sample_info[1, 7]
  if (ncol(sample_info) == 8)  reference_genome <- sample_info[1, 8]
  if (ncol(sample_info) == 9){ reference_genome <- sample_info[1, 8]; simple_repeat_list <- sample_info[1, 9] }
} else {
  adapter_2 <- sample_info[1, 6]
  organism  <- sample_info[1, 7]
  panel     <- sample_info[1, 8]
  if (ncol(sample_info) == 9)  reference_genome <- sample_info[1, 9]
  if (ncol(sample_info) == 10){ reference_genome <- sample_info[1, 9]; simple_repeat_list <- sample_info[1,10] }
}
if (is.na(adapter_2) || adapter_2 == "") adapter_2 <- adapter_1

# Normalise paths that are used by external tools
mutation_file <- normalizePath(mutation_file, mustWork = TRUE)
bam_file_in   <- normalizePath(bam_file_in, mustWork = TRUE)
if (exists("reference_genome")) {
  reference_genome <- normalizePath(reference_genome, mustWork = TRUE)
}

# ----------------- Paths for intermediates -----------------

tmp_dir <- file.path(wd, "tmp")
if (!dir.exists(tmp_dir)) dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

bam_file_slim     <- file.path(tmp_dir, paste0(sample_name, ".SLIM.bam"))
bam_file_slim_bai <- paste0(bam_file_slim, ".bai")
bed_path          <- file.path(tmp_dir, paste0(sample_name, ".regions.bed"))

# ----------------- Genome & chromosome naming -----------------

ref_genome <- fun_load_genome(organism)
chr_no     <- fun_load_chr_no(organism)

if (ref_genome@user_seqnames[[1]] == "chr1") {
  chromosomes <- paste0("chr", c(seq_len(chr_no - 2), "X", "Y"))
} else {
  chromosomes <- paste0("", c(seq_len(chr_no - 2), "X", "Y"))
}

# ----------------- Load mutations -----------------

df_mutation <- fun_load_mutation(mutation_file, sample_name, ref_genome, chr_no)

# Normalize required columns & fallbacks
req_cols <- c("Sample","Mut_type","Chr","Pos","Ref","Alt","SimpleRepeat_TRF","Neighborhood_sequence")
missing <- setdiff(req_cols, colnames(df_mutation))
if (length(missing) > 0) {
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

# ----------------- Build BED around mutations (±200; merge if next within 400bp) -----------------

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

# ----------------- samtools pipeline (streamed) -----------------

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
    qpath(tmp_dir), threads, qpath(bam_file_in), qpath(bed_path),
    threads, mem_per_thread, qpath(sort_prefix), qpath(bam_file_slim)
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
    qpath(tmp_dir), threads, qpath(bam_file_in), qpath(reference_genome), qpath(bed_path),
    threads, mem_per_thread, qpath(sort_prefix), qpath(bam_file_slim)
  )
} else {
  stop(sprintf("Unsupported input extension: %s", input_ext))
}

message("Trimming and sorting reads into slim BAM (tmp/)...")
run_cmd(cmd)
run_cmd(sprintf("bash -lc 'set -euo pipefail; samtools index %s'", qpath(bam_file_slim)))

# ----------------- load slim BAM -----------------

df_bam <- fun_load_bam(bam_file_slim)

# ----------------- Filter mutations to BAM contigs -----------------

targets <- names(Rsamtools::scanBamHeader(bam_file_slim)[[1]]$targets)
missing_chr <- setdiff(unique(df_mutation$Chr), targets)
if (length(missing_chr) > 0) {
  message(sprintf("[%s] Skipping contigs not in BAM: %s",
                  sample_name, paste(missing_chr, collapse = ", ")))
}
df_mutation <- dplyr::filter(df_mutation, Chr %in% targets)

# If nothing left after subsetting, emit empty result and exit
if (nrow(df_mutation) == 0) {
  message(sprintf("[%s] No mutations left on BAM contigs after subsetting; writing empty result.", sample_name))
  con <- gzfile(out_path, "w")
  writeLines("", con)  # empty gz file
  close(con)
  message(sprintf("MicroSEC finished successfully (empty result). Output: %s", out_path))
  q(save = "no", status = 0)
}

# ----------------- MicroSEC core (single sample) -----------------

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

msec            <- result[[1]]
homology_search <- result[[2]]
mut_depth       <- result[[3]]

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

# ----------------- Save final result -----------------

fun_save(msec, out_path)
message(sprintf("MicroSEC finished successfully. Output: %s", out_path))
