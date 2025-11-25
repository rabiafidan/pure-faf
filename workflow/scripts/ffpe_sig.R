##############################################
## 0. LOAD PACKAGES & REFERENCE GENOME
##############################################

library(tidyverse)
library(MutationalPatterns)
library(pheatmap)

## Choose BSgenome based on Snakemake param
if (snakemake@params[["ref_genome"]] == "grch38") {
  library(BSgenome.Hsapiens.UCSC.hg38)
  ref_genome <- BSgenome.Hsapiens.UCSC.hg38      # BSgenome object for mut_matrix()
  genome_name <- "BSgenome.Hsapiens.UCSC.hg38"   # string for read_vcfs_as_granges()
} else if (snakemake@params[["ref_genome"]] == "grch37") {
  library(BSgenome.Hsapiens.UCSC.hg19)
  ref_genome <- BSgenome.Hsapiens.UCSC.hg19
  genome_name <- "BSgenome.Hsapiens.UCSC.hg19"
}


##############################################
## 1. INPUTS FROM SNAKEMAKE
##############################################

## Each of these is a *vector* of VCF paths (one per sample)
unfiltered_vcf                <- snakemake@input[["unfiltered_vcf"]]
qual_filtered_vcf             <- snakemake@input[["qual_filtered_vcf"]]
vault_filtered_vcf            <- snakemake@input[["vault_filtered_vcf"]]
vault_plus_sr_filtered_vcf    <- snakemake@input[["vault_plus_sr_filtered_vcf"]]
msec_all_filtered_vcf         <- snakemake@input[["msec_all_filtered_vcf"]]
VAF_filtered_vault_vcf        <- snakemake@input[["VAF_filtered_vault_vcf"]]
VAF_filtered_vault_plus_sr_vcf<- snakemake@input[["VAF_filtered_vault_plus_sr_vcf"]]
VAF_filtered_msec_all_vcf     <- snakemake@input[["VAF_filtered_msec_all_vcf"]]

## Sample names are shared across all filtering steps
sample_names <- snakemake@params[["sample_names"]]
ref_gen      <- snakemake@params[["ref_genome"]]   # string ("grch37"/"grch38")

##############################################
## 2. LOAD & FORMAT FFPE SIGNATURES
##############################################
ffpe_sig <- read.csv("resources/ffpesig.csv")

## Convert to MutationalPatterns-style mutation type: e.g.
## "ACG>AT" -> "A[C>A]T"
mut_type <- paste0(
  substr(ffpe_sig$MutationType, 5, 5), "[",  # flanking base right of central
  substr(ffpe_sig$MutationType, 1, 3), "]",  # central triplet with > (e.g. C>A)
  substr(ffpe_sig$MutationType, 7, 7)        # flanking base left of central
)

## For safety, store this in the data.frame
ffpe_sig$MP_type <- mut_type

##############################################
## 3. BUILD A NAMED LIST OF VCF VECTORS BY STEP
##############################################

## Each element is a character vector of VCF paths (one per sample)
vcf_by_step <- list(
  unfiltered                 = unfiltered_vcf,
  qual_filtered              = qual_filtered_vcf,
  vault_filtered             = vault_filtered_vcf,
  vault_plus_sr_filtered     = vault_plus_sr_filtered_vcf,
  msec_all_filtered          = msec_all_filtered_vcf,
  VAF_filtered_vault         = VAF_filtered_vault_vcf,
  VAF_filtered_vault_plus_sr = VAF_filtered_vault_plus_sr_vcf,
  VAF_filtered_msec_all      = VAF_filtered_msec_all_vcf
)

filter_steps <- names(vcf_by_step)

##############################################
## 4. HELPER: READ VCFs & MAKE 96-CONTEXT MATRIX
##############################################

make_mut_matrix_for_step <- function(vcf_paths, sample_names, genome_name, ref_genome) {
  # Read VCFs as GRanges, one GRanges per sample
  vcfs_as_gr <- read_vcfs_as_granges(
    vcf_files   = vcf_paths,
    sample_names = sample_names,
    genome       = genome_name
  )
  
  # Build 96-trinucleotide mutation count matrix
  # rows: 96 mutation contexts, cols: samples
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

## We now align ffpe_sig rows to the *exact* order of mut_types_ref
idx <- match(mut_types_ref, ffpe_sig$MP_type)

## If any mutation context from MutationalPatterns is missing in FFPE sigs,
## stop and complain
if (any(is.na(idx))) {
  missing_types <- mut_types_ref[is.na(idx)]
  stop("Some mutation types from mut_matrix are missing in ffpe_sig:\n",
       paste(missing_types, collapse = ", "))
}

ffpe_sig_aligned <- ffpe_sig[idx, ]
rownames(ffpe_sig_aligned) <- mut_types_ref

## Select the numeric signature columns.
## Adjust this if your CSV has different column names.
## Here I assume the 2nd and 3rd columns are the signatures.
ffpe_mat <- as.matrix(ffpe_sig_aligned[, 2:3])

## Give more informative names to the signatures
## (adjust to match your actual column names / meaning)
colnames(ffpe_mat) <- c("FFPE_Repaired", "FFPE_Unrepaired")

## Normalise each signature to sum to 1 (probability distribution over 96 contexts)
ffpe_mat_norm <- sweep(ffpe_mat, 2, colSums(ffpe_mat), FUN = "/")

##############################################
## 7. FIT FFPE SIGNATURES TO EACH STEP (CONTRIBUTIONS)
##############################################

## For each filtering step, we fit the FFPE signatures to the mutation matrix
## using non-negative least squares via fit_to_signatures().
##
## fit_to_signatures() returns:
## - contribution:  (#signatures x #samples) matrix of counts
## - reconstructed: 96 x samples matrix reconstructed from signatures
##
## We store the `contribution` matrices in a list.

contrib_list <- lapply(mut_mats, function(mat_step) {
  fit <- fit_to_signatures(
    mut_matrix = mat_step,
    signatures = ffpe_mat_norm
  )
  # fit$contribution has rows = signatures, cols = samples
  fit$contribution
})

##############################################
## 8. BUILD SUMMARY MATRICES (ROWS = STEPS, COLS = SAMPLES)
##    FOR EACH FFPE SIGNATURE
##############################################

## We'll create two matrices:
## - contrib_unrepaired[step, sample]
## - contrib_repaired[step, sample]
##
## These will be filled from contrib_list[[step]][signature, sample]

n_steps   <- length(filter_steps)
n_samples <- length(sample_names)

contrib_unrepaired <- matrix(
  NA_real_,
  nrow = n_steps,
  ncol = n_samples,
  dimnames = list(filter_steps, sample_names)
)

contrib_repaired <- contrib_unrepaired

for (step_name in filter_steps) {
  contrib_step <- contrib_list[[step_name]]  # rows = signatures, cols = samples
  
  ## Add contributions to the appropriate row
  contrib_unrepaired[step_name, ] <- contrib_step["FFPE_Unrepaired", ]
  contrib_repaired[step_name, ]   <- contrib_step["FFPE_Repaired", ]
}

##############################################
## 9. OPTIONALLY CONVERT TO FRACTION OF MUTATIONS
##############################################

## If you prefer fractions instead of raw counts, we can divide by the total
## number of mutations per sample *per step*.

contrib_unrepaired_frac <- contrib_unrepaired
contrib_repaired_frac   <- contrib_repaired

for (step_name in filter_steps) {
  mat_step <- mut_mats[[step_name]]
  total_muts_per_sample <- colSums(mat_step)  # total mutations in that step+sample
  
  contrib_unrepaired_frac[step_name, ] <- contrib_unrepaired[step_name, ] /
    total_muts_per_sample
  
  contrib_repaired_frac[step_name, ] <- contrib_repaired[step_name, ] /
    total_muts_per_sample
}

##############################################
## 10. HEATMAPS (ROWS = FILTERING STEPS, COLS = SAMPLES)
##############################################

## Heatmap of *fraction* of mutations attributed to FFPE_Unrepaired
h1=pheatmap(
  contrib_unrepaired_frac,
  main           = "FFPE Unrepaired signature - fraction of mutations contributing",
  cluster_rows   = FALSE,  # keep steps in the specified order
  cluster_cols   = FALSE,  # keep samples in the specified order
  display_numbers = TRUE,  # show values in the cells
  number_format  = "%.2f",
  color = colorRampPalette(c("steelblue", "white", "indianred4"))(100),
  breaks = seq(0, 1, length.out = 101)
)

## Heatmap of *fraction* of mutations attributed to FFPE_Repaired
h2=pheatmap(
  contrib_repaired_frac,
  main           = "FFPE Repaired signature - fraction of mutations contributing",
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  display_numbers = TRUE,
  number_format  = "%.2f",
  color = colorRampPalette(c("steelblue", "white", "indianred4"))(100),
  breaks = seq(0, 1, length.out = 101)
)

##############################################
## 11.SAVE
##############################################
pdf(snakemake@output[["plot"]], width=10, height=4)
h1
h2
dev.off()

write.table(contrib_repaired_frac, file=snakemake@output[["repaired_table"]], sep="\t", quote=FALSE)
write.table(contrib_unrepaired_frac, file=snakemake@output[["unrepaired_table"]], sep="\t", quote=FALSE)