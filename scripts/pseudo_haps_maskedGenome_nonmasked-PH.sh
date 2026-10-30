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

# maskedGenome_masked-PH ->  Mapping: Masked Genome; PH: Masked Genome
# maskedGenome_nonmasked-PH -> Mapping: Masked Genome; PH: Regualar Genome
# maskedGenome-noRealign_masked-PH -> Mapping: Masked no realign Genome; PH: Masked Genome
# maskedGenome-noRealign_nonmasked-PH -> Mapping: Masked no realign Genome; PH: Regualar Genome

TAG=maskedGenome_nonmasked-PH

# which mapping run to use:
#   $BASE/PopGLen             masked genome, with realignment
#   $BASE/PopGLen_norealign   masked genome, no realignment
POPGLEN_DIR=/home/rbrennan/rw_edna_wgs/PopGLen
 
# whichref for pseudohaplotypes
#   GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.SNPs.masked.fasta     masked
#   GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.fasta                 non-masked
REF_FASTA=/home/rbrennan/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.fasta
 
####################################################################

BASE=/home/rbrennan/rw_edna_wgs
OUTDIR=$BASE/analysis/pseudo_haps

BAM_LIST="$OUTDIR/bam_list_${TAG}.txt"
OUT_PREFIX="$OUTDIR/low_cov_maf05_${TAG}_pseudohap"

ls "$POPGLEN_DIR"/results/mapping/bams/*b_ricei_masked.rmdup.realn.clip.bam > "$BAM_LIST"

#samtools faidx /home/rbrennan/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.SNPs.masked.fasta

# change chr ids to match:
#awk 'NR==FNR{map[$1]=$2; next} {if ($1 in map) $1=map[$1]; print}' OFS='\t' \
#  ~/rw_edna_wgs/scripts/chrom_map.tsv ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.txt > ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.renamed.txt

#angsd sites index ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.renamed.txt

angsd -bam "$BAM_LIST" \
  -ref  "$REF_FASTA"\
  -sites ~/rw_edna_wgs/analysis/high_cov/high_cov_sites_maf05.renamed.txt \
  -rf ~/rw_edna_wgs/analysis/pseudo_haps/chrs.txt \
  -doHaploCall 1 \
  -doCounts 1 \
  -minMapQ 30 \
  -minQ 30 \
  -remove_bads 1 \
  -uniqueOnly 1 \
  -nThreads 10 \
  -minInd 0 \
  -setMinDepthInd 1 \
  -out "$OUT_PREFIX"

zcat "${OUT_PREFIX}.haplo.gz" > "${OUT_PREFIX}.haplo"

awk 'NR==FNR{map[$2]=$1; next} {if ($1 in map) $1=map[$1]; print}' OFS='\t'  \
  ~/rw_edna_wgs/scripts/chrom_map.tsv "${OUT_PREFIX}.haplo" \
  > "${OUT_PREFIX}.renamed.haplo"
