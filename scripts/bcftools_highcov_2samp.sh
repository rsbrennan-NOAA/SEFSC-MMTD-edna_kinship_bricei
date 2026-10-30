#!/bin/bash
#SBATCH --job-name=bam_geno
#SBATCH --time=3-00:00:00
#SBATCH --mem=44G
#SBATCH -c 4
#SBATCH -n 1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard
 
module load bio/bcftools/1.21 
module load bio/htslib/1.19

# break if errors
set -euo pipefail

REF=~/reference_genomes/rices_whale/GCF_028023285.1_mBalRic1.hap2_genomic.fna
VCF=~/rw_edna_wgs/analysis/high_cov/bricei_only_variants_chr.maf05.vcf.gz
CHROM_MAP=~/rw_edna_wgs/scripts/chrom_map.tsv
BAMDIR=~/rw_edna_wgs/PopGLen/results/mapping/bams
OUTDIR=~/rw_edna_wgs/analysis/variants/
# ---------------

cd "$OUTDIR"

#sites file from the VCF, with chr names converted to NW_ names (reference/BAM)
bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' "$VCF" |
  awk 'NR==FNR{map[$1]=$2; next} {if ($1 in map) $1=map[$1]; print}' OFS='\t' \
    "$CHROM_MAP" - |
  bgzip -c > sites.tsv.gz
tabix -f -s1 -b2 -e2 sites.tsv.gz

# break if any chr names were not converted
if zcat sites.tsv.gz | cut -f1 | grep -qv '^NW_'; then
  echo "ERROR: some chromosome names were not converted to NW_ names" >&2
  exit 1
fi

# call variants at those sites, same alleles, keep hom calls
bcftools mpileup -f "$REF" -T sites.tsv.gz -a AD,DP \
    "$BAMDIR"/Bric026.b_ricei.rmdup.realn.clip.bam \
    "$BAMDIR"/Bric032.b_ricei.rmdup.realn.clip.bam \
    "$BAMDIR"/Bric045.b_ricei.rmdup.realn.clip.bam \
    "$BAMDIR"/Bric046.b_ricei.rmdup.realn.clip.bam \
    -Ou |
  bcftools call -m -C alleles -T sites.tsv.gz --threads 4 -Oz -o bam_indivs.vcf.gz

bcftools index bam_indivs.vcf.gz
