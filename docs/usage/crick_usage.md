# Usage from within the Francis Crick Institute

**NOTE: you need to run the main snakemake job on the login node. You need to use tmux or screen sessions for this purpose.** 
**In order to reattach a session, you need to be on the same login node where you created the session (login006 or login007).**

You can use the environment module system to load singularity (and conda if not already loaded):

### 1. load anaconda:
```bash
ml load Anaconda3/2024.10-1
conda init
```

### 2. you can then create and activate a conda environment for snakemake:

```bash
conda config --set channel_priority strict
conda create --name pure-faf -c conda-forge -c bioconda snakemake cyvcf2 snakemake-executor-plugin-slurm

# or you can pip install snakemake-executor-plugin-slurm after creating the conda environment:
# pip install snakemake-executor-plugin-slurm
```

### 3. create a tmux or screen session and activate the conda environment:
```bash
screen -S pure-faf_session
conda activate pure-faf
```
### 4. load singularity:
```bash
ml load Singularity/3.11.3
```

Then follow the instructions in the [general usage guide](general_usage.md) to clone the repository, configure and run the workflow.
To see an example of snakemake profile configuration for slurm, see [turajlic lab usage](turajlic_lab_usage.md).
