# Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/rabiafidan/pure-faf).


If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository or its DOI.

### 0) Dependencies
You need `snakemake`, `cyvcf2`, `conda` and `singularity` to run this workflow.

```bash
conda create --name pure-faf -c conda-forge -c bioconda snakemake singularity cyvcf2
conda activate pure-faf
```
If you are running this workflow on an HPC environment, make sure to also install the correct executor plugin. For example:

```
pip install snakemake-executor-plugin-slurm
```
### 1) clone the repository
```bash
git clone https://github.com/rabiafidan/pure-faf.git
```
To run the workflow from command line, change the working directory.

```bash
cd path/to/pure-faf
```

### 2) configure the pipeline
Adjust options in the default config file `config/config.yaml`.


### 3) run the pipeline
Before running the complete workflow, you can perform a dry run using:

```bash
snakemake --dry-run --use-conda --use-singularity
```
To run the workflow locally with test files using **conda** and **singularity**:

```bash
snakemake --cores 2 --use-conda --use-singularity
```
you can add `--no-temp` flag if you need to keep every intermidate file (for debugging purposes).

You can also run the workflow on an **HPC environment** using an executor plugin. See snakemake documentation for more details on how to set up executor plugins.
See [turajlic lab usage](turajlic_lab_usage.md) for an example of slurm profile configuration.