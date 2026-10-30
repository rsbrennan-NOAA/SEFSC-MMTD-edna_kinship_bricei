library(tidyverse)

# ibis analysis

ibd1 <- read.csv("analysis/ibis/bricei_only_variants_chr.maf05_ibd1.coef", header=T, sep="\t")
ibd2 <- read.csv("analysis/ibis/bricei_only_variants_chr.maf05_ibd2.coef", header=T, sep="\t")

colnames(ibd1) <- c("sample1","sample2","kinship_coefficient","IBD2_fraction","segment_count","degree_of_relatedness")
colnames(ibd2) <- c("sample1","sample2","kinship_coefficient","IBD2_fraction","segment_count","degree_of_relatedness")

ibd_seg1 <- read.csv("analysis/ibis/bricei_only_variants_chr.maf05_ibd1.seg", header=T, sep="\t")
ibd_seg2 <- read.csv("analysis/ibis/bricei_only_variants_chr.maf05_ibd2.seg", header=T, sep="\t")

colnames(ibd_seg1) <- c("sample1","sample2","chrom","phys_start_pos","phys_end_pos","IBD_type","genetic_start_pos","genetic_end_pos","genetic_seg_length","marker_count","error_count","error_density")
colnames(ibd_seg2) <- c("sample1","sample2","chrom","phys_start_pos","phys_end_pos","IBD_type","genetic_start_pos","genetic_end_pos","genetic_seg_length","marker_count","error_count","error_density")

p_coef <- ggplot(ibd1, aes(kinship_coefficient)) +
              geom_histogram(alpha = 1, bins=60)  +
  theme_bw()

p_coef2 <- ggplot(ibd2, aes(kinship_coefficient)) +
              geom_histogram(alpha = 1, bins=60)  +
  theme_bw()

ibd1[which(ibd1$kinship_coefficient > 0.24),]

# in macaq paper, plot n IBD1 > 12cM by sum IBD1 > 12cM
# To estimate realized relatedness (hereafter called rIBD), we summed the total length of shared IBD segments ≥ 8 cM per dyad (sumIBD1 + 2*sumIBD2) and divided it by the total diploid length of the macaque genome (2,970 cM) to calculate the percentage of shared DNA between individuals.
ibd_seg1$genetic_seg_length

length_threshold_cM <- 7      
total_genome_cM     <- 2773   

rIBD_table <- ibd_seg2 |>
  filter(genetic_seg_length >= length_threshold_cM) |>
  mutate(dyad = paste(pmin(sample1, sample2), pmax(sample1, sample2), sep = "_")) |>
  group_by(dyad) |>
  summarise(
    sumIBD1 = sum(genetic_seg_length[IBD_type == "IBD1"]),
    sumIBD2 = sum(genetic_seg_length[IBD_type == "IBD2"]),
    .groups = "drop"
  ) |>
  mutate(r_IBD = (sumIBD1 + 2 * sumIBD2) / total_genome_cM)

rIBD_table

ggplot(rIBD_table, aes(r_IBD)) +
              geom_histogram(alpha = 1, bins=60)  +
  theme_bw()


# 1. Genome coverage fraction, not just total length. A true parent-offspring pair (or duplicate/self-comparison) should have IBD segments spanning essentially the entire genome with almost no gaps — that's a categorically different, saturating pattern from "generally elevated background with real recombination breaks visible throughout." Background inflation makes segments longer/more numerous everywhere, but it's much harder for it to fake near-100% genome coverage in one specific pair.

# 2. IBD2 fraction. Parent-offspring dyads should have ~0% IBD2 (you can never inherit both copies from one parent). Full siblings should sit around 25%. This ratio is a relative signature within a pair, not an absolute scale — it survives inflation better than total r_IBD does.

# 3. Segment count relative to coverage. For a given amount of genome covered, more numerous/shorter segments suggest a more distant true relationship whose long apparent segments got merged by low marker density, vs. fewer/longer segments suggesting a genuinely close, recently-descended relationship.

dyad_features <- ibd_seg2 |>
  filter(genetic_seg_length >= length_threshold_cM) |>   # <-- restored; this was missing
  mutate(dyad = paste(pmin(sample1, sample2), pmax(sample1, sample2), sep = "_")) |>
  group_by(dyad) |>
  summarise(
    sumIBD1 = sum(genetic_seg_length[IBD_type == "IBD1"]),
    sumIBD2 = sum(genetic_seg_length[IBD_type == "IBD2"]),
    n_segments = n(),
    max_seg = max(genetic_seg_length),
    .groups = "drop"
  ) |>
  mutate(
    coverage_fraction = (sumIBD1 + sumIBD2) / total_genome_cM,
    ibd2_fraction = sumIBD2 / (sumIBD1 + sumIBD2),
    r_IBD = (sumIBD1 + 2 * sumIBD2) / total_genome_cM
  ) |>
  arrange(desc(r_IBD))

dyad_features

as.data.frame(dyad_features)

# find the max length segment for each.
dyad_features$r_IBD

fullRel <- read.csv("analysis/pseudo_haps/relatednessNewData.csv", header=T)

fullRel$Indv1 <- gsub("Bricei", "Bric", fullRel$Indv1)
fullRel$Indv2 <- gsub("Bricei", "Bric", fullRel$Indv2)

dyad_features <- dyad_features|>
  separate(dyad, into = c("Indv1", "Indv2"), sep = "_", remove = FALSE)
# use rab

dyad_features <- dyad_features |>
  mutate(id1 = pmin(Indv1, Indv2),
         id2 = pmax(Indv1, Indv2))

fullRel <- fullRel |>
  mutate(id1 = pmin(Indv1, Indv2),
         id2 = pmax(Indv1, Indv2))

merged <- dyad_features |>
  left_join(fullRel, by = c("id1", "id2"), suffix = c("", "_full"))

ggplot(data=merged, aes(rab, r_IBD)) +
  geom_point(size=2) +
  theme_bw()

merged <- merged |>
  arrange(desc(r_IBD))
merged[c("dyad", "max_seg", "r_IBD", "rab")]
tail(merged[c("dyad", "max_seg", "r_IBD", "rab")])

ibd_seg1 |>
  ggplot(aes(x = markers_per_cM)) +
  geom_histogram(bins = 60) +
  scale_x_log10()

ibd_seg1 <- ibd_seg1 |>
  mutate(markers_per_cM = marker_count / genetic_seg_length)

ibd_seg1 |>
  ggplot(aes(x = markers_per_cM)) +
  geom_histogram(bins = 60) +
  scale_x_log10()

density_threshold <- 80

ibd_seg1_filtered <- ibd_seg1 |>
  filter(markers_per_cM >= density_threshold)

dyad_features_filtered <- ibd_seg1_filtered |>
  filter(genetic_seg_length >= length_threshold_cM) |>
  mutate(dyad = paste(pmin(sample1, sample2), pmax(sample1, sample2), sep = "_")) |>
  group_by(dyad) |>
  summarise(
    sumIBD1 = sum(genetic_seg_length[IBD_type == "IBD1"]),
    sumIBD2 = sum(genetic_seg_length[IBD_type == "IBD2"]),
    n_segments = n(),
    max_seg = max(genetic_seg_length),
    .groups = "drop"
  ) |>
  mutate(
    ibd2_fraction = sumIBD2 / (sumIBD1 + sumIBD2),
    r_IBD = (sumIBD1 + 2 * sumIBD2) / total_genome_cM
  ) |>
  arrange(desc(r_IBD))




ibd_seg1_filtered |>
  mutate(markers_per_cM = marker_count / genetic_seg_length) |>
  ggplot(aes(x = genetic_seg_length, y = markers_per_cM)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  scale_x_log10() +
  ylim(0,620)





merged_filtered <- dyad_features_filtered |>
  left_join(select(merged, dyad, rab), by = "dyad")

cor(merged_filtered$r_IBD, merged_filtered$rab, use = "complete.obs")
cor(merged$r_IBD, merged$rab, use = "complete.obs")   # your original, for comparison


ggplot(data=merged_filtered, aes(rab, r_IBD)) +
  geom_point(size=2) +
  theme_bw()

merged <- merged |>
  arrange(desc(r_IBD))
merged[c("dyad", "max_seg", "r_IBD", "rab")]
tail(merged[c("dyad", "max_seg", "r_IBD", "rab")])


ibd_seg1 <- ibd_seg1 |>
  mutate(markers_per_cM = marker_count / genetic_seg_length)

ibd_seg1 |>
  ggplot(aes(x = markers_per_cM)) +
  geom_histogram(bins = 60) +
  scale_x_log10()


merged_filtered <- merged_filtered |>
  mutate(markers_per_cM = marker_count / genetic_seg_length)


merged_filtered |>
  ggplot(aes(x = rab, y = r_IBD)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "r_IBD vs rab (filtered)")

merged_filtered |>
  ggplot(aes(x = markers_per_cM)) +
  geom_histogram(bins = 60) +
  scale_x_log10()

dyad_marker_density <- ibd_seg1_filtered |>
  mutate(dyad = paste(pmin(sample1, sample2), pmax(sample1, sample2), sep = "_")) |>
  group_by(dyad) |>
  summarise(
    mean_markers_per_cM = mean(markers_per_cM),
    median_markers_per_cM = median(markers_per_cM),
    .groups = "drop"
  )

merged_filtered <- merged_filtered |>
  left_join(dyad_marker_density, by = "dyad")

merged_filtered |>
  ggplot(aes(x = markers_per_cM)) +
  geom_histogram(bins = 60) +
  scale_x_log10()



ibd_seg1 |>
  mutate(markers_per_cM = marker_count / genetic_seg_length) |>
  ggplot(aes(x = genetic_seg_length, y = markers_per_cM)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess") +
  scale_x_log10()
















relatedness_highcov
which.min( relatedness_highcov$pmr)

relatedness_highcov <- relatedness_highcov |>
  separate(pair, into = c("Indv1", "Indv2"), sep = " - ", remove = FALSE)


relatedness_highcov[1:4,1:4]


fullRel <- fullRel |>
  mutate(id1 = pmin(Indv1, Indv2),
         id2 = pmax(Indv1, Indv2))

relatedness_highcov <- relatedness_highcov |>
  mutate(id1 = pmin(Indv1, Indv2),
         id2 = pmax(Indv1, Indv2))

merged <- relatedness_highcov |>
  left_join(fullRel, by = c("id1", "id2"), suffix = c("", "_full"))

#pmin()/pmax() alphabetically sort each row's pair so Bric001-Bric002 and Bric002-Bric001 both become id1 = Bric001, id2 = Bric002, letting the join match regardless of original order.

#Note: this will produce duplicate Indv1/Indv2 columns from both data frames (with _full suffix added to fullRel's version
library(ggpubr)
rab <- ggplot(data=merged, aes(rab, pmr)) +
  geom_point(size=2) +
  theme_bw()

merged$rab
as.data.frame(merged[24,])

merged$pmr[is.na(merged$rab)]


R1_R0 <- ggplot(merged, aes(x=R1, y=R0)) +
              geom_point(alpha=0.5, size=2)

R1_R0 <- ggplot(merged, aes(x = R1, y = R0, color = pmr < 0.24)) +
              geom_point(alpha = 0.5, size = 2)

merged <- merged |>
  mutate(group = case_when(
    pmr < 0.25 & rab > 0.3  ~ "Both",
    pmr < 0.25 & rab <= 0.3 ~ "pmr < 0.24 only",
    pmr >= 0.25 & rab > 0.3 ~ "rab > 0.4 only",
    TRUE                     ~ "Neither"
  ))|>
  arrange(desc(group == "Neither")) 

R1_R0 <- ggplot(merged, aes(x = R1, y = R0, color = group)) +
              geom_point(alpha = 1, size = 2) +
              scale_color_manual(values = c("Both" = "purple",
                                             "pmr < 0.24 only" = "orange",
                                             "rab > 0.4 only" = "blue",
                                             "Neither" = "grey70")) +
  theme_bw()
R1_R0
ggarrange(rab,R1_R0)