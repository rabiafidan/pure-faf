log <- file(snakemake@log[[1]], open = "at")
sink(log, append = TRUE)
sink(log, type = "message", append = TRUE)
message("-----------------------------------------------")
message("-----------------------------------------------")
message(date())

##############################################
## 0. INSTALL AND LOAD PACKAGES & REFERENCE GENOME
##############################################
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

options(repos = BiocManager::repositories())

if (!require("MutationalPatterns", quietly = TRUE)) {
  BiocManager::install("MutationalPatterns")
}

# pheatmap is on CRAN
if (!require("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}

library(tidyverse)
library(MutationalPatterns)
library(pheatmap)

## Choose BSgenome based on Snakemake param
if (snakemake@params[["ref_genome"]] == "grch38") {
  if (!require("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
  }
  library(BSgenome.Hsapiens.UCSC.hg38)
  ref_genome <- BSgenome.Hsapiens.UCSC.hg38      # BSgenome object for mut_matrix()
  genome_name <- "BSgenome.Hsapiens.UCSC.hg38"   # string for read_vcfs_as_granges()
} else if (snakemake@params[["ref_genome"]] == "grch37") {
  if (!require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  }
  library(BSgenome.Hsapiens.UCSC.hg19)
  ref_genome <- BSgenome.Hsapiens.UCSC.hg19
  genome_name <- "BSgenome.Hsapiens.UCSC.hg19"
} else {
  stop("ref_genome must be 'grch37' or 'grch38'")
}

##############################################
## 1. INPUTS FROM SNAKEMAKE
##############################################

## Each of these is a *vector* of VCF paths (one per sample)
unfiltered_vcf                 <- snakemake@input[["unfiltered_vcf"]]
qual_filtered_vcf              <- snakemake@input[["qual_filtered_vcf"]]
vault_filtered_vcf             <- snakemake@input[["vault_filtered_vcf"]]
vault_plus_sr_filtered_vcf     <- snakemake@input[["vault_plus_SR_filtered_vcf"]]
msec_all_filtered_vcf          <- snakemake@input[["msec_all_filtered_vcf"]]
VAF_filtered_vault_vcf         <- snakemake@input[["VAF_filtered_vault_vcf"]]
VAF_filtered_vault_plus_sr_vcf <- snakemake@input[["VAF_filtered_vault_plus_SR_vcf"]]
VAF_filtered_msec_all_vcf      <- snakemake@input[["VAF_filtered_msec_all_vcf"]]

## Sample names are shared across all filtering steps
sample_names <- snakemake@params[["sample_names"]]
ref_gen      <- snakemake@params[["ref_genome"]]   # string ("grch37"/"grch38")

##############################################
## 2. LOAD & FORMAT FFPE SIGNATURES
##############################################
## FFPE signature CSV is assumed to have:
## - a column 'MutationType' like "ACG>AT"
## - numeric columns for the two FFPE signatures (here assumed columns 2 and 3)

ffpe_sig <- read.csv("resources/ffpesig.csv", stringsAsFactors = FALSE)

## Convert to MutationalPatterns-style mutation type: e.g. "ACG>AT" -> "A[C>A]T"
mut_type <- paste0(
  substr(ffpe_sig$MutationType, 5, 5), "[",  # base right of central
  substr(ffpe_sig$MutationType, 1, 3), "]",  # central triplet with >
  substr(ffpe_sig$MutationType, 7, 7)        # base left of central
)

ffpe_sig$MP_type <- mut_type

##############################################
## 3. BUILD A NAMED LIST OF VCF VECTORS BY STEP
##############################################
## Helper: from many VAF-filtered paths -> list(step_name -> vector of paths)
build_vaf_step_list <- function(paths, sample_names, prefix_label) {
  df <- tibble(
    path = paths,
    base = basename(paths)
  ) %>%
    mutate(
      tum = sub(
        paste0("^2-(.+?)_", prefix_label, "_VAF_.*$"),
        "\\1",
        base
      ),
      thr = sub(
        "^.*_VAF_([0-9.]+)_filtered\\.vcf\\.gz$",
        "\\1",
        base
      )
    )
  
  ## Hard fail if parsing breaks
  if (any(is.na(df$tum)) || any(is.na(df$thr))) {
    stop("Failed to parse tumour names or VAF thresholds for ", prefix_label)
  }
  
  steps <- list()
  
  for (thr in sort(unique(df$thr))) {
    df_thr <- df[df$thr == thr, ]
    
    ## Strict sample check
    if (!setequal(df_thr$tum, sample_names)) {
      stop(
        "Sample mismatch for VAF ", thr, " in ", prefix_label, "\n",
        "Expected: ", paste(sample_names, collapse = ", "), "\n",
        "Found: ", paste(df_thr$tum, collapse = ", ")
      )
    }
    
    ## Reorder paths to match sample_names EXACTLY
    vec <- df_thr$path[match(sample_names, df_thr$tum)]
    
    step_name <- paste0(prefix_label, "_VAF_", thr)
    steps[[step_name]] <- vec
  }
  
  steps
}

## Base (non-VAF) steps
base_steps <- list(
  unfiltered             = unfiltered_vcf,
  qual_filtered          = qual_filtered_vcf,
  vault_filtered         = vault_filtered_vcf,
  vault_plus_SR_filtered = vault_plus_sr_filtered_vcf,
  msec_all_filtered      = msec_all_filtered_vcf
)

## Build per-VAF steps for each filter type
vault_vaf_steps <- build_vaf_step_list(
  paths        = VAF_filtered_vault_vcf,
  sample_names = sample_names,
  prefix_label = "vault_filtered"
)

vault_sr_vaf_steps <- build_vaf_step_list(
  paths        = VAF_filtered_vault_plus_sr_vcf,
  sample_names = sample_names,
  prefix_label = "vault_plus_SR_filtered"
)

msec_vaf_steps <- build_vaf_step_list(
  paths        = VAF_filtered_msec_all_vcf,
  sample_names = sample_names,
  prefix_label = "msec_all_filtered"
)

## Combine everything: names become rows in the heatmaps
vcf_by_step <- c(
  base_steps,
  vault_vaf_steps,
  vault_sr_vaf_steps,
  msec_vaf_steps
)
for (nm in names(vcf_by_step)) {
  message(nm, ":\n  ", paste(basename(vcf_by_step[[nm]]), collapse = "\n  "))
}

filter_steps <- names(vcf_by_step)

##############################################
## 4. HELPER: READ VCFs & MAKE 96-CONTEXT MATRIX
##############################################

make_mut_matrix_for_step <- function(vcf_paths, sample_names, genome_name, ref_genome) {
  # Read VCFs as GRanges, one GRanges per sample
  vcfs_as_gr <- read_vcfs_as_granges(
    vcf_files    = vcf_paths,
    sample_names = sample_names,
    genome       = genome_name
  )
  
  # Build 96-trinucleotide mutation count matrix
  mut_matrix(
    vcf_list   = vcfs_as_gr,
    ref_genome = ref_genome
  )
}

##############################################
## 5. BUILD MUTATION MATRICES FOR ALL STEPS
##############################################

mut_mats <- lapply(vcf_by_step, make_mut_matrix_for_step,
                   sample_names = sample_names,
                   genome_name  = genome_name,
                   ref_genome   = ref_genome)

## Use the first matrix as reference for row order (mutation types)
mut_types_ref <- rownames(mut_mats[[1]])

## Check that all matrices have identical rownames
stopifnot(all(vapply(mut_mats,
                     function(m) identical(rownames(m), mut_types_ref),
                     logical(1))))

##############################################
## 6. ALIGN FFPE SIGNATURES TO MUTATIONALPATTERNS ROW ORDER
##############################################

idx <- match(mut_types_ref, ffpe_sig$MP_type)

if (any(is.na(idx))) {
  missing_types <- mut_types_ref[is.na(idx)]
  stop("Some mutation types from mut_matrix are missing in ffpe_sig:\n",
       paste(missing_types, collapse = ", "))
}

ffpe_sig_aligned <- ffpe_sig[idx, ]
rownames(ffpe_sig_aligned) <- mut_types_ref

## Extract the two FFPE signatures (assume columns 2 and 3 are numeric)
ffpe_mat <- as.matrix(ffpe_sig_aligned[, 2:3])

## Give informative names (adjust if needed to match your CSV)
colnames(ffpe_mat) <- c("FFPE_Repaired", "FFPE_Unrepaired")

## Normalise each FFPE signature to sum to 1 over 96 contexts
ffpe_mat_norm <- sweep(ffpe_mat, 2, colSums(ffpe_mat), FUN = "/")

##############################################
## 7. COSINE SIMILARITY: SAMPLE PROFILES vs FFPE SIGNATURES
##############################################
## For each filtering step:
##  - build sample profiles (96-context fractions per sample)
##  - compute cosine similarity between FFPE signatures and sample profiles
##
## Result: two matrices
##  - sim_unrepaired[step, sample]
##  - sim_repaired[step, sample]
##############################################

n_steps   <- length(filter_steps)
n_samples <- length(sample_names)

sim_unrepaired <- matrix(
  NA_real_,
  nrow = n_steps,
  ncol = n_samples,
  dimnames = list(filter_steps, sample_names)
)

sim_repaired <- sim_unrepaired

## store variant counts per step x sample
variant_counts <- matrix(
  NA_integer_,
  nrow = n_steps,
  ncol = n_samples,
  dimnames = list(filter_steps, sample_names)
)

for (step_name in filter_steps) {
  mat_step <- mut_mats[[step_name]]  # 96 x samples
  
  # Total mutations per sample, to normalise columns
  total_muts_per_sample <- colSums(mat_step)
  
  # store counts for later use in the plot
  variant_counts[step_name, ] <- total_muts_per_sample
  
  # Avoid division by zero: for samples with zero mutations, keep NA
  zero_cols <- total_muts_per_sample == 0
  
  # Build sample profiles: each column = fraction of mutations in each context
  sample_profiles <- sweep(mat_step, 2, total_muts_per_sample, FUN = "/")
  
  # Set columns with zero mutations to NA (no defined profile)
  if (any(zero_cols)) {
    sample_profiles[, zero_cols] <- NA_real_
  }
  
  ## Compute cosine similarity between FFPE signatures and sample profiles.
  ## Both have rows = mutation types, columns = entity (signature/sample).
  ## cos_sim_matrix() returns a matrix with:
  ##   rows = colnames(first argument)
  ##   cols = colnames(second argument)
  sim_mat <- cos_sim_matrix(ffpe_mat_norm, sample_profiles)
  # sim_mat["FFPE_Repaired", sample] is the similarity we want
  
  sim_unrepaired[step_name, ] <- sim_mat["FFPE_Unrepaired", ]
  sim_repaired[step_name, ]   <- sim_mat["FFPE_Repaired", ]
}

## build label matrices "cosine - #variants"
numbers_unrepaired <- matrix(
  NA_character_,
  nrow = n_steps,
  ncol = n_samples,
  dimnames = dimnames(sim_unrepaired)
)

numbers_repaired <- numbers_unrepaired

for (i in seq_len(n_steps)) {
  for (j in seq_len(n_samples)) {
    cs_u <- sim_unrepaired[i, j]
    cs_r <- sim_repaired[i, j]
    nvar <- variant_counts[i, j]
    
    if (!is.na(cs_u) && !is.na(nvar)) {
      numbers_unrepaired[i, j] <- sprintf("%.2f (%d snvs)", cs_u, as.integer(nvar))
    }
    if (!is.na(cs_r) && !is.na(nvar)) {
      numbers_repaired[i, j] <- sprintf("%.2f (%d snvs)", cs_r, as.integer(nvar))
    }
  }
}

##############################################
## 8. HEATMAPS (ROWS = FILTERING STEPS, COLS = SAMPLES)
##    VALUES = COSINE SIMILARITY (0..1)
##############################################

h1 <- pheatmap(
  sim_unrepaired,
  main            = "Cosine similarity to FFPE Unrepaired signature\n(labels: cosine similarity (number of snvs))",
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = numbers_unrepaired,
  number_color    = "black",
  color           = colorRampPalette(c("steelblue", "white", "indianred4"))(100),
  breaks          = seq(0, 1, length.out = 101)
)

h2 <- pheatmap(
  sim_repaired,
  main            = "Cosine similarity to FFPE Repaired signature\n(labels: cosine similarity (number of snvs))",
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  display_numbers = numbers_repaired,
  number_color    = "black",
  color           = colorRampPalette(c("steelblue", "white", "indianred4"))(100),
  breaks          = seq(0, 1, length.out = 101)
)

##############################################
## 9. SAVE
##############################################
pdf(snakemake@output[["plot1"]], width = 10, height = 4)
h1
dev.off()

pdf(snakemake@output[["plot2"]], width = 10, height = 4)
h2
dev.off()

## Save the similarity matrices as tables
write.table(
  sim_repaired,
  file  = snakemake@output[["repaired_table"]],
  sep   = "\t",
  quote = FALSE
)

write.table(
  sim_unrepaired,
  file  = snakemake@output[["unrepaired_table"]],
  sep   = "\t",
  quote = FALSE
)

# Close sinks and log file connection
sink(type = "message")
sink()
close(log)