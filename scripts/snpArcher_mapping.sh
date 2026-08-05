#!/bin/bash
#SBATCH --job-name=snparcher_controller
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard


source ~/.bashrc
conda activate snparcher

snakemake \
  --executor slurm \
  --snakefile ~/snpArcher/workflow/Snakefile \
  --directory ~/rw_edna_wgs/snpArcher \
  --workflow-profile ~/rw_edna_wgs/snpArcher/workflow-profile \
  --use-conda \
  --conda-prefix ~/snparcher_envs \
  --jobs 15 \
  --retries 3 \
  --latency-wait 20 \
  --until combine_qc_metrics mosdepth 
