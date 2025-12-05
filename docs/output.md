# Output files and structure

The main output of the workflow are the formalin filtered VCF files located in `results/odir/filtered_vcf/`.

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