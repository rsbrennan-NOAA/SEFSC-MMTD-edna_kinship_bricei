#!/bin/bash

#SBATCH --job-name=pseudohap   ## job name
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH --cpus-per-task 10   
#SBATCH --mem=90G
#SBATCH -t 16:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 
#SBATCH --ntasks=1 
#SBATCH --partition standard

##########################################################################################
module load bio/angsd/0.940
module load bio/samtools/1.23

BASE=/home/rbrennan/rw_edna_wgs
OUTDIR=$BASE/analysis/pseudo_haps

ls "$BASE"/PopGLen/results/mapping/bams/*b_ricei.rmdup.realn.clip.bam > "$OUTDIR/bam_list.txt"

#samtools faidx /home/rbrennan/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.fasta

angsd sites index /home/rbrennan/rw_edna_wgs/analysis/pseudo_haps/high_cov_sites.renamed.txt

angsd -bam "$OUTDIR/bam_list.txt" \
  -ref "/home/rbrennan/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.fasta" \
  -sites "$OUTDIR/high_cov_sites.renamed.txt" \
  -doMajorMinor 3 \
  -doHaploCall 1 \
  -doCounts 1 \
  -minMapQ 30 \
  -minQ 30 \
  -remove_bads 1 \
  -uniqueOnly 1 \
  -nThreads 10 \
  -minInd 0 \
  -out "$OUTDIR/low_cov_pseudohap"
