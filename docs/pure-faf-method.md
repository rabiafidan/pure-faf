<img width="2560" height="1440" alt="pure-faf_workflow" src="https://github.com/user-attachments/assets/c470b105-4ee2-48bd-a638-5fe70578a8bb" />

## Workflow overview

This workflow is a best-practice workflow for variant filtering of formalin fixed tumour samples.
The workflow is built using [snakemake](https://snakemake.readthedocs.io/en/stable/) and consists of the following steps:

1. Annotating rhe VCF with a custom panel of normals (PON)
This PON has been created from RepSeq samples. Beta distribution testing is performed to identify variants that likely to come from PON.
2. Applying a general quality filter with creteria including PON, AD and ROQ
3. Splitting variants into confidence tiers using an AD threshold
4. Intersecting low-confidence indels with a second tool's indel calls. Intersection is indel type-specific (insertion or deletion), and only position-based (alt allele doesn't matter).
5. Running microsec twice with different parameters for low-confidence and high-confidence tiers
6. Merging the results back annotating VCF with microsec results
7. Applyying different filtering methods on the annotated VCF
8. Applying a post-filter VAF filter 
9. Quality control using FFPEsig repaired and unrepaired signatures

