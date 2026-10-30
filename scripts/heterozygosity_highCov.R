library(snpR)
library(tidyverse)

# Stage 1: read VCF into a snpRdata object
dat <- read_vcf("analysis/high_cov/bricei_only_variants_chr.maf05.vcf.gz")

# Get sample names, then re-import with an ID column in the sample metadata
ids <- colnames(genotypes(dat))

dat <- read_vcf("analysis/high_cov/bricei_only_variants_chr.maf05.vcf.gz", sample.meta = data.frame(ind = ids))

dat <- calc_ho(dat, facets = "ind")

ind_ho <- get.snpR.stats(dat, facets = "ind", stats = "ho")
ind_ho$weighted.means

#-------------------------
library(vcfR)
library(tidyverse)

library(vcfR)
library(tidyverse)

vcf <- read.vcfR("analysis/high_cov/bricei_only_variants_chr.maf05.vcf.gz")
vcf <- vcf[is.biallelic(vcf), ]

gt <- extract.gt(vcf, element = "GT")
ad <- extract.gt(vcf, element = "AD")

ab_df <- tibble(
  sample = rep(colnames(gt), each = nrow(gt)),
  gt = as.vector(gt),
  ad = as.vector(ad)
) |>
  filter(!is.na(gt), !is.na(ad), !str_detect(ad, "\\.")) |>
  separate_wider_delim(ad, delim = ",", names = c("ref", "alt")) |>
  mutate(across(c(ref, alt), as.numeric),
         depth = ref + alt,
         ab = alt / depth) |>
  filter(depth >= 15)

ab_df |>
  filter(gt %in% c("0/1", "0|1", "1|0")) |>
  ggplot(aes(ab)) +
  geom_histogram(bins = 50) +
  facet_wrap(~sample, scales = "free_y")



#--------------------------------------------

vcf_bam2 <- read.vcfR("~/Documents/projects/rw_edna_wgs/analysis/variants/bam_indivs.vcf.gz")
vcf_bam2 <- vcf_bam[is.biallelic(vcf_bam), ]

gt <- extract.gt(vcf_bam2, element = "GT")
ad <- extract.gt(vcf_bam2, element = "AD")

ab_bam <- tibble(
  sample = rep(colnames(gt), each = nrow(gt)),
  gt = as.vector(gt),
  ad = as.vector(ad)
) |>
  filter(!is.na(gt), !is.na(ad), !str_detect(ad, "\\.")) |>
  separate_wider_delim(ad, delim = ",", names = c("ref", "alt"),
                       too_few = "align_start") |>
  mutate(across(c(ref, alt), as.numeric),
         alt = replace_na(alt, 0),
         depth = ref + alt,
         ab = alt / depth) |>
  filter(depth >= 15)

ab_bam |>
  filter(gt %in% c("0/1", "0|1", "1|0")) |>
  ggplot(aes(ab)) +
  geom_histogram(bins = 70) +
  facet_wrap(~sample, scales = "free_y")



ab_bam |>
  filter(sample == "Bric026",
         gt %in% c("0/1", "0|1", "1|0")) |>
  ggplot(aes(ab)) +
  geom_histogram(bins = 70) +
  labs(title = "Bric026", x = "alt allele fraction", y = "count")




ab_bam |>
  filter(sample == "Bric026",
         gt %in% c("0/0", "0|0", "1/1", "1|1")) |>
  mutate(minor = pmin(ab, 1 - ab)) |>
  summarise(
    n_hom = n(),
    mean_depth = mean(depth),
    prop_with_minor = mean(minor > 0),
    mean_minor = mean(minor)
  )




hom_summary <- function(df) {
  df |>
    filter(gt %in% c("0/0", "0|0", "1/1", "1|1")) |>
    mutate(minor = pmin(ab, 1 - ab)) |>
    group_by(sample) |>
    summarise(
      n_hom = n(),
      mean_depth = mean(depth),
      prop_with_minor = mean(minor > 0),
      mean_minor = mean(minor)
    )
}

hom_compare <- bind_rows(
  hom_summary(ab_df) |> mutate(source = "Original VCF"),
  hom_summary(ab_bam |> filter(sample == "Bric026")) |> mutate(source = "New Sample")
) |>
  arrange(desc(prop_with_minor))

as.data.frame(hom_compare)

