library(tidyverse)

compute_ref_alt_freq <- function(geno_path, ind_path, group_label) {
  ind_ids <- read_table(ind_path, col_names = FALSE)$X1
  geno_lines <- readLines(geno_path)
  geno_mat <- do.call(rbind, strsplit(geno_lines, ""))  # rows = SNPs, cols = individuals
  colnames(geno_mat) <- ind_ids

  # geno: 2 = ref, 0 = alt, 9 = missing
  tibble(sample_id = ind_ids, group = group_label) |>
    rowwise() |>
    mutate(
      n_ref = sum(geno_mat[, sample_id] == "2"),
      n_alt = sum(geno_mat[, sample_id] == "0"),
      n_called = n_ref + n_alt,
      frac_ref = n_ref / n_called,
      frac_alt = n_alt / n_called
    ) |>
    ungroup()
}

highcov_freqs <- compute_ref_alt_freq(
  "analysis/pseudo_haps/highcov_27samples_merged.geno",
  "analysis/pseudo_haps/highcov_27samples_merged.ind",
  "HighCov"
)

lowcov_freqs <- compute_ref_alt_freq(
  "analysis/pseudo_haps/lowcov_filtered.geno",
  "analysis/pseudo_haps/lowcov_filtered.ind",
  "LowCov"
)

all_freqs <- bind_rows(highcov_freqs, lowcov_freqs)

highcov_mean_ref <- mean(highcov_freqs$frac_ref, na.rm = TRUE)

p <- ggplot(all_freqs, aes(x = sample_id, y = frac_ref, fill = group)) +
  geom_col() +
  geom_hline(yintercept = highcov_mean_ref, linetype = "dashed", color = "firebrick3") +
  annotate("text", x = 1, y = highcov_mean_ref + 0.03,
           label = paste0("high-cov mean = ", round(highcov_mean_ref, 2)),
           hjust = 0, size = 3, color = "firebrick3") +
  labs(y = "Proportion of ref allele calls", x = NULL, 
       title = "Reference bias") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p

ggsave("figures/reference_bias_maskedGenome_masked-PH.png", p, w=7, h=4)

# -----------------------
# 
# 
# -----------------
# loop over all combos:



library(tidyverse)

compute_ref_alt_freq <- function(geno_path, ind_path, group_label) {
  ind_ids <- read_table(ind_path, col_names = FALSE)$X1
  geno_lines <- readLines(geno_path)
  geno_mat <- do.call(rbind, strsplit(geno_lines, ""))  # rows = SNPs, cols = individuals
  colnames(geno_mat) <- ind_ids

  # geno: 2 = ref, 0 = alt, 9 = missing
  tibble(sample_id = ind_ids, group = group_label) |>
    rowwise() |>
    mutate(
      n_ref = sum(geno_mat[, sample_id] == "2"),
      n_alt = sum(geno_mat[, sample_id] == "0"),
      n_called = n_ref + n_alt,
      frac_ref = n_ref / n_called,
      frac_alt = n_alt / n_called
    ) |>
    ungroup()
}


TAGS <- c(
  "nonmaskedGenome_nonmasked-PH",
  "maskedGenome_masked-PH",
  "maskedGenome_nonmasked-PH",
  "maskedGenome-noRealign_masked-PH",
  "maskedGenome-noRealign_nonmasked-PH"

)


highcov_freqs_by_tag <- list()
lowcov_freqs_by_tag <- list()

for (tag in TAGS) {
  highcov_freqs <- compute_ref_alt_freq(
    paste0("analysis/pseudo_haps/highcov_", tag, "_27samples_merged.geno"),
    paste0("analysis/pseudo_haps/highcov_", tag, "_27samples_merged.ind"),
    "HighCov"
  )

  lowcov_freqs <- compute_ref_alt_freq(
    paste0("analysis/pseudo_haps/lowcov_filtered_", tag, ".geno"),
    paste0("analysis/pseudo_haps/lowcov_filtered_", tag, ".ind"),
    "LowCov"
  )

  cat(sprintf("[%s] high-cov samples: %d, low-cov samples: %d\n",
              tag, nrow(highcov_freqs), nrow(lowcov_freqs)))

  highcov_freqs_by_tag[[tag]] <- highcov_freqs
  lowcov_freqs_by_tag[[tag]] <- lowcov_freqs
}

head(lowcov_freqs_by_tag[["maskedGenome_nonmasked-PH"]])
head(lowcov_freqs_by_tag[["maskedGenome-noRealign_nonmasked-PH"]])



highcov_freqs_by_tag <- list()
lowcov_freqs_by_tag <- list()

for (tag in TAGS) {
  highcov_freqs <- compute_ref_alt_freq(
    paste0("analysis/pseudo_haps/highcov_", tag, "_27samples_merged.geno"),
    paste0("analysis/pseudo_haps/highcov_", tag, "_27samples_merged.ind"),
    "HighCov"
  )

  lowcov_freqs <- compute_ref_alt_freq(
    paste0("analysis/pseudo_haps/lowcov_filtered_", tag, ".geno"),
    paste0("analysis/pseudo_haps/lowcov_filtered_", tag, ".ind"),
    "LowCov"
  )

  cat(sprintf("[%s] high-cov samples: %d, low-cov samples: %d\n",
              tag, nrow(highcov_freqs), nrow(lowcov_freqs)))

  highcov_freqs_by_tag[[tag]] <- highcov_freqs
  lowcov_freqs_by_tag[[tag]] <- lowcov_freqs
}




# ---- Combine high-cov + low-cov, compute high-cov mean ref frequency, per combo ----

all_freqs_by_tag <- list()
highcov_mean_ref_by_tag <- list()

for (tag in TAGS) {
  highcov_freqs <- highcov_freqs_by_tag[[tag]]
  lowcov_freqs <- lowcov_freqs_by_tag[[tag]]

  all_freqs <- bind_rows(highcov_freqs, lowcov_freqs)

  highcov_mean_ref <- mean(highcov_freqs$frac_ref, na.rm = TRUE)

  cat(sprintf("[%s] high-cov mean ref frac: %.3f\n", tag, highcov_mean_ref))

  all_freqs_by_tag[[tag]] <- all_freqs
  highcov_mean_ref_by_tag[[tag]] <- highcov_mean_ref
}


TAG_LABELS <- list(
  "nonmaskedGenome_nonmasked-PH" = "Non-Masked genome; PH w/ non-masked",
  "maskedGenome_masked-PH" = "Masked genome; PH w/ masked",
  "maskedGenome_nonmasked-PH" = "Masked genome; PH w/ non-masked",
  "maskedGenome-noRealign_masked-PH" = "Masked no-realign genome; PH w/ masked",
  "maskedGenome-noRealign_nonmasked-PH" = "Masked no-realign genome; PH w/ non-masked"
)

# plot

plots_by_tag <- list()

for (tag in TAGS) {
  all_freqs <- all_freqs_by_tag[[tag]]
  highcov_mean_ref <- highcov_mean_ref_by_tag[[tag]]

  p <- ggplot(all_freqs, aes(x = sample_id, y = frac_ref, fill = group)) +
    geom_col() +
    geom_hline(yintercept = highcov_mean_ref, linetype = "dashed", color = "firebrick3") +
    annotate("text", x = 1, y = highcov_mean_ref + 0.03,
             label = paste0("high-cov mean = ", round(highcov_mean_ref, 2)),
             hjust = 0, size = 3, color = "firebrick3") +
    labs(y = "Proportion of ref allele calls", x = NULL,
         title = paste0("Reference bias: ", TAG_LABELS[[tag]])) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  print(p)

  ggsave(paste0("figures/reference_bias_", tag, ".png"), p, w = 7, h = 4)

  plots_by_tag[[tag]] <- p
}


ref_bias_plots_by_tag <- list()
non_missing_plots_by_tag <- list()

for (tag in TAGS) {
  all_freqs <- all_freqs_by_tag[[tag]]
  lowcov_freqs <- lowcov_freqs_by_tag[[tag]]
  highcov_mean_ref <- highcov_mean_ref_by_tag[[tag]]

  # ---- reference bias plot (unchanged) ----
  p_bias <- ggplot(all_freqs, aes(x = sample_id, y = frac_ref, fill = group)) +
    geom_col() +
    geom_hline(yintercept = highcov_mean_ref, linetype = "dashed", color = "firebrick3") +
    annotate("text", x = 1, y = highcov_mean_ref + 0.03,
             label = paste0("high-cov mean = ", round(highcov_mean_ref, 2)),
             hjust = 0, size = 3, color = "firebrick3") +
    labs(y = "Proportion of ref allele calls", x = NULL,
         title = paste0("Reference bias: ", TAG_LABELS[[tag]])) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  print(p_bias)
  ggsave(paste0("figures/reference_bias_", tag, ".png"), p_bias, w = 7, h = 4)
  ref_bias_plots_by_tag[[tag]] <- p_bias

  # ---- non-missing sites plot  ----
  p_missing <- ggplot(lowcov_freqs, aes(x = sample_id, y = n_called)) +
    geom_col(fill = "steelblue") +
    labs(y = "Non-missing sites", x = NULL,
         title = paste0("Non-missing sites: ", TAG_LABELS[[tag]])) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ylim(0,10000)

  print(p_missing)
  ggsave(paste0("figures/non_missing_sites_", tag, ".png"), p_missing, w = 7, h = 4)
  non_missing_plots_by_tag[[tag]] <- p_missing
}






