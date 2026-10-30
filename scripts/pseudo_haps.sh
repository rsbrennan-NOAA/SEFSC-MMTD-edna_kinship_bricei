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

ls "$BASE"/PopGLen/results/mapping/bams/*b_ricei_masked.rmdup.realn.clip.bam > "$OUTDIR/bam_list.txt"

#samtools faidx /home/rbrennan/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.SNPs.masked.fasta

# change chr ids to match:
awk 'NR==FNR{map[$1]=$2; next} {if ($1 in map) $1=map[$1]; print}' OFS='\t' \
  ~/rw_edna_wgs/scripts/chrom_map.tsv ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.txt |\
 sort -k1  > ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.renamed.txt

angsd sites index ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.renamed.txt

angsd -bam "$OUTDIR/bam_list.txt" \
  -ref "/home/rbrennan/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.fasta" \
  -sites ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.renamed.txt \
  -doHaploCall 1 \
  -doCounts 1 \
  -minMapQ 30 \
  -minQ 20 \
  -remove_bads 1 \
  -uniqueOnly 1 \
  -nThreads 10 \
  -minInd 0 \
  -setMinDepthInd 1 \
  -out "$OUTDIR/low_cov_maf05_pseudohap_noMajMin"

zcat ~/rw_edna_wgs/analysis/pseudo_haps/low_cov_maf05_pseudohap_noMajMin.haplo.gz > ~/rw_edna_wgs/analysis/pseudo_haps/low_cov_maf05_pseudohap_noMajMin.haplo

awk 'NR==FNR{map[$2]=$1; next} {if ($1 in map) $1=map[$1]; print}' OFS='\t'  \
  ~/rw_edna_wgs/scripts/chrom_map.tsv ~/rw_edna_wgs/analysis/pseudo_haps/low_cov_maf05_pseudohap_noMajMin.haplo  \
   > ~/rw_edna_wgs/analysis/pseudo_haps/low_cov_maf05_pseudohap_noMajMin.renamed.haplo
