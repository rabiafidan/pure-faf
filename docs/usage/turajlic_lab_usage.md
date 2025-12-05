# Usage from within the Turajlic Lab

**NOTE: you need to run the main snakemake job on the login node. You need to use tmux or screen sessions for this purpose.** 
**In order to reattach a session, you need to be on the same login node where you created the session (login006 or login007).**

You can use the environment module system to load singularity (and conda if not already loaded):

1. load anaconda:
```bash
ml load Anaconda3/2024.10-1
conda init
```

2. you can then create and activate a conda environment for snakemake:

```bash
conda config --set channel_priority strict
conda create --name pure-faf -c conda-forge -c bioconda snakemake cyvcf2 snakemake-executor-plugin-slurm

# or you can pip install snakemake-executor-plugin-slurm after creating the conda environment:
# pip install snakemake-executor-plugin-slurm
```
3. clone the repository
```bash
git clone https://github.com/rabiafidan/pure-faf.git
```
To run the workflow from command line, change the working directory.

```bash
cd path/to/pure-faf
```

4. configure the pipeline
Adjust options in the default config file `config/config.yaml`.

5. create a tmux or screen session and activate the conda environment:
```bash
screen -S pure-faf_session
conda activate pure-faf
```

6. load singularity:
```bash
ml load Singularity/3.11.3
```
7. run the pipeline

Before running the complete workflow, you can perform a dry run using:
```bash
snakemake --dry-run --executor slurm --profile /nemo/project/proj-tracerX/working/PIPELINES/Snakemake -p
```
To run the workflow
```bash
snakemake --executor slurm --profile /nemo/project/proj-tracerX/working/PIPELINES/Snakemake -p
```

you can add `--no-temp` flag if you need to keep every intermidate file (for debugging purposes).

### profile configuration
Our Turajlic lab profile configuration is as pasted below for reference. For other labs, you need to change your slurm_account name and you need to change the singularity-args to make sure it binds all the necessary paths for your lab.

```bash
cat /nemo/project/proj-tracerX/working/PIPELINES/Snakemake/config.yaml 
```
```yaml
executor: slurm
latency-wait: 80
rerun-triggers: mtime
keep-going: true
rerun-incomplete: true
singularity-args: '--bind /tmp,/flask/reference/Genomics/aws-igenomes/Homo_sapiens/,/camp/project/proj-tracerX/,/nemo/project/proj-tracerX/,/nemo/ess-manifest/proj-turajlics-rcc,/nemo/project/proj-turajlics-wgs --env GITHUB_PAT=${GITHUB_PAT}'
jobs: 70
slurm-keep-successful-logs: true
slurm-delete-logfiles-older-than: 0
slurm-requeue: true
slurm-efficiency-report: true
slurm-efficiency-report-path: results/.efficiency_reports/
use-conda: true
use-singularity: true  
use-envmodules: true

default-resources:
    slurm_partition: "ncpu"
    slurm_account: "u_turajlics"
    runtime: 5000
    mem_mb_per_cpu: 4000
    slurm_cpus_per_task: 1 
```
