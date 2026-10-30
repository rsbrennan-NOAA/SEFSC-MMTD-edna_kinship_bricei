#!/bin/bash
#SBATCH --job-name=mask
#SBATCH --time=03:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard


module load bio/bedtools/2.31.1 

cd ~/rw_edna_wgs/analysis/high_cov

# remanme chr ids to match
echo "renaming chromosome ids"

zcat bricei_only_variants_chr.maf01.vcf.gz | \
awk 'NR==FNR{map[$1]=$2; next} /^#/{print; next} {if ($1 in map) $1=map[$1]; print}' \
  OFS='\t' ~/rw_edna_wgs/scripts/chrom_map.tsv - > bricei_only_variants_chr.maf01.renamed.vcf

#convert to  BED to avoid memory crash
echo "convert vcf to bed"
awk 'BEGIN{OFS="\t"} !/^#/{print $1, $2-1, $2}' bricei_only_variants_chr.maf01.renamed.vcf > variants.bed

cd /home/rbrennan/reference_genomes/rices_whale/

# mask reference using the mitogenome version
echo "mask reference"
bedtools maskfasta \
  -fi /home/rbrennan/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.fasta \
  -bed ~/rw_edna_wgs/analysis/high_cov/variants.bed \
  -fo GCF_028023285.1_mBalRic1.hap2_genomic_mitogenome.SNPs.masked.fasta

