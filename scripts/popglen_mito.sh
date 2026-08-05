#!/bin/bash
#SBATCH --job-name=popglen_mito
#SBATCH --time=7-00:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard
 
 
source ~/.bashrc
conda activate popglen
 
snakemake \
  --executor slurm \
  --snakefile ~/rw_edna_wgs/PopGLen_mito/workflow/Snakefile \
  --directory ~/rw_edna_wgs/PopGLen_mito \
  --workflow-profile ~/rw_edna_wgs/PopGLen_mito/profiles/sedna \
  --use-conda \
  --conda-prefix ~/popglen_envs \
  --use-singularity \
  --singularity-args "--bind /home/rbrennan/reference_genomes" \
  --singularity-prefix ~/popglen_singularity \
  --jobs 20 \
  --retries 3 \
  --latency-wait 20 \
  --rerun-incomplete \
  --until bam_clipoverlap 
