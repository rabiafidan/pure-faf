# Input data
### You need these following files:
1. Mutect2 VCF files (with AD, DP and ROQ annotations)
2. Strelka2 indel VCF (or any other secondary variant caller should theoratically work)
3. Sample bam/cram file
4. Human reference genome fasta you used for alignment and variant calling

### You need the following information:
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
