#!/bin/bash
#SBATCH --job-name=sites_hicov
#SBATCH --time=05:00:00
#SBATCH --mem=24G
#SBATCH -c 4
#SBATCH -n 1
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o /home/rbrennan/rw_edna_wgs/logs/%x_%j.out
#SBATCH -e /home/rbrennan/rw_edna_wgs/logs/%x_%j.err
#SBATCH --partition standard

module load bio/bcftools
module load bio/angsd/0.940

cd ~/rw_edna_wgs/analysis/high_cov
# grab sample names, filter to those starting with "Bric"
bcftools query -l briceiOutgroupAll_PASS_variants_chr.vcf.gz | grep '^Bric' > bric_samples.txt

bcftools view -S bric_samples.txt -O z -o bricei_only_variants_chr.vcf.gz briceiOutgroupAll_PASS_variants_chr.vcf.gz

# refilter vcf:

bcftools +fill-tags bricei_only_variants_chr.vcf.gz -- -t AC,AN,AF | \
  bcftools view -m2 -M2 -v snps -q 0.01:minor -O z -o bricei_only_variants_chr.maf01.vcf.gz

tabix -p vcf bricei_only_variants_chr.maf01.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' bricei_only_variants_chr.maf01.vcf.gz > high_cov_sites.txt

angsd sites index high_cov_sites.txt

# and get maf 0.05- might want for downstream pmr

bcftools +fill-tags bricei_only_variants_chr.vcf.gz -- -t AC,AN,AF | \
  bcftools view -m2 -M2 -v snps -q 0.05:minor -O z -o bricei_only_variants_chr.maf05.vcf.gz

tabix -p vcf bricei_only_variants_chr.maf05.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' bricei_only_variants_chr.maf05.vcf.gz > high_cov_sites_maf05.txt


