## Workflow overview

This workflow is a best-practice workflow for variant filtering of formalin fixed tumour samples.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Annotating rhe VCF with a custom panel of normals (PON)
2. Applying a general quality filter with creteria including PON, AD and ROQ
3. Splitting variants into confidence tiers using an AD threshold
4. Intersecting low-confidence indels with a second tool's indel calls. Intersection is indel type-specific (insertion or deletion), and only position-based (alt allele doesn't matter).
5. Running microsec twice with different parameters for low-confidence and high-confidence tiers
6. Merging the results back annotating VCF with microsec results
7. Applyying different filtering methods on the annotated VCF
8. Applying a post-filter VAF filter 
9. Quality control using FFPEsig repaired and unrepaired signatures

See the [worflow diagram](README.md).

## Running the workflow

### Input data
You need these following files:
1. Mutect2 VCF files (with AD, DP and ROQ annotations)
2. Strelka2 indel VCF (or any other secondary variant caller should theoratically work)
3. Sample bam/cram file
4. Human reference genome fasta you used for alignment and variant calling

You need the following information:
1. tumour and normal sample names
2. Sequencing read length
3. Sequencing adapters
 
Using the input files and information, you need to create a sample sheet as follows:

|tumour_name          |normal_name          |Mutect2_vcf                                                                                    |Strelka2_indel_vcf                                                                                                                                                                |tumour_bam_cram                                                                                                             |sequencing_read_length|sequencing_adapter_1|sequencing_adapter_2|
|---------------------|---------------------|-----------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------|----------------------|--------------------|--------------------|
|K1058_K1058_HB005_D02|K1058_K1058_BC001_D01|path/to/ROQ_flagged_K1058_HB005_D02.vcf.gz                                                     |path/to/K1058_HB005_D02_vs_K1058_BC001_D01/K1058_HB005_D02_vs_K1058_BC001_D01.strelka.somatic_indels.vcf.gz                                                                       |path/to/K1058_HB005_D02/K1058_HB005_D02.recal.cram                                                                          |300                   |AGCGAGAT-CCAGTATC   |AGCGAGAT-CCAGTATC   |
|K1058_K1058_HB003_D02|K1058_K1058_BC001_D01|path/to/ROQ_flagged_K1058_HB003_D02.vcf.gz                                                     |path/to/K1058_HB003_D02_vs_K1058_BC001_D01/K1058_HB003_D02_vs_K1058_BC001_D01.strelka.somatic_indels.vcf.gz                                                                       |path/to/K1058_HB003_D02/K1058_HB003_D02.recal.cram                                                                          |300                   |CTTGCTAG-TCATCTCC   |CTTGCTAG-TCATCTCC   |



### Parameters

You then need to configure the workflow by updating the [config/config.yaml](config/config.yaml).
This table lists all parameters that can be used to run the workflow.
All parameters except exclude_samples are mandatory.

| parameter                         | type | details                                                                                           | VAULT setup                    |
| --------------------------------- | ---- | ------------------------------------------------------------------------------------------------- | ------------------------------ |
| **samplesheet**                   |      |                                                                                                   |                                |
| path                              | str  | path to samplesheet.                                                                              |                                |
| **ref_genome**                    |      |                                                                                                   |                                |
| ref_genome_version                | str  | choose one of GRCh38 or GRCh37, case-insensitive.                                                 |                                |
| ref_genome_fasta                  | str  | path to the human reference genome.                                                               |                                |
| **VCF quality filter thresholds** |      |                                                                                                   |                                |
| AD                                | num  | minimum alternative allele count                                                                  | 4                              |
| ROQ                               | num  | ROQ threshold for read orientation quality                                                        | 20                             |
| DP                                | num  | minimum read depth                                                                                | 0 - no seperate DP filtering   |
| **PON parameters**                |      |                                                                                                   |                                |
| PON_alpha                         | num  | Beta test significance level                                                                      | 0.05                           |
| **Microsec parameters**           |      |                                                                                                   |                                |
| AD_tier_threshold                 | num  | maximum AD for low-confidence tier                                                                | 10                             |
| pval_threshold_high_AD            | num  | pvalue threshold for `fun_analysis` function for high-confidence tier                             | 1e-4                           |
| pval_threshold_low_AD             | num  | pvalue threshold for `fun_analysis` function for low-confidence tier                              | 1e-6                           |
| **Post-filter thresholds**        |      |                                                                                                   |                                |
| VAF                               | array| a list of VAF values you want to apply after formalin filtering. All will appear in the QC plot.  |                                |
| **Other parameters**              |      |                                                                                                   |                                |
| odir                              | str  | output directory. This will appear under results/                                                 |                                |
| exclude_samples                   | array| a list of sample names you want to exclude from the analysis. Default is empty list.              |                                |                   

### Output
```bash
results/odir/
    ├── logs
        ├── 1-step1_sample1.log
        ├── 2-step2_sample1.log
        ├── ...
    ├── filtered_vcf #MAIN OUTPUT FOLDER
        ├── 1-sample1_msec_all_filtered.vcf.gz #passing all microsec filters
        ├── 1-sample2_vault_filtered.vcf.gz #filtering creteria used in VAULT paper
        ├── 1-sample1__vault_plus_SR_filtered.vcf.gz #VAULT createria + simple repeat filter (relevant if whole genome data, otherwise almost no difference)
        ├── 2-..... #various VAF filtered VCFs
        ├── 2-..... #various VAF filtered VCFs
    ├── microsec_input # sample info and mutation info files for microsec
    ├── microsec_output #microsec tsv output files for samples and tiers
    ├── ffpe_sig #QC plots and tables using FFPEsig repaired and unrepaired signatures
    ├── temp # intermediate files. Most is automatically deleted after workflow completion. Remaining ones may be useful or can be deleted manually.
    └── PON # custom RepSeq panel of normals annotation and beta distribution testing results. Kind of intermediate, can be deleted if not needed.
```