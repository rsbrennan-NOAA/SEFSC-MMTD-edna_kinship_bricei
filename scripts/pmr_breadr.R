library(tidyverse)
library(BREADR)


# analyze the haps from the low cov samples


low_haps <- read.table("analysis/pseudo_haps/low_cov_pseudohap.haplo",
  header=T)

ind_summary <- low_haps |>
  select(starts_with("ind")) |>
  summarise(across(everything(), ~ sum(. != "N"))) |>
  pivot_longer(everything(), names_to = "individual", values_to = "n_called") |>
  mutate(pct_called = 100 * (n_called / nrow(low_haps)))

ind_summary


processEigenstratNoFilter <- function (indfile, genofile, snpfile, filter_length = NULL, pop_pattern = NULL, 
    filter_deam = FALSE, outfile = NULL, chromosomes = NULL, 
    verbose = TRUE, byPop = FALSE, minMerge = 2) 
{
    if (!file.exists(indfile)) {
        stop(paste0("indfile ", indfile, " does not exist!"))
    }
    if (!file.exists(genofile)) {
        stop(paste0("genofile ", genofile, " does not exist!"))
    }
    if (!file.exists(snpfile)) {
        stop(paste0("snpfile ", snpfile, " does not exist!"))
    }
    if ((!is.numeric(minMerge)) | (minMerge < 2)) {
        stop("minMerge must be a number >= 2.")
    }
    if (is.null(filter_length)) {
        filter_length <- 1e+05
        if (verbose) {
            cat(sprintf("No site distance filter.\nUsing default minimum of 1e5.\n"))
        }
    }
    if (!is.numeric(filter_length) | (filter_length <= 0)) {
        stop("filter_length must be a positive number.")
    }
    ind_raw <- read_ind(indfile)
    if (!all(pop_pattern %in% ind_raw$pop)) {
        paste0("Following populations in pop_pattern do not exist in indfile: ", 
            paste0(setdiff(pop_pattern, ind_raw$pop), collapse = ", ")) %>% 
            stop()
    }
    if (!is.logical(filter_deam)) {
        stop("filter_deam must be a logical (boolean) variable, i.e., TRUE or FALSE.")
    }
    if (!is.null(outfile)) {
        if (!file.exists(dirname(outfile))) {
            stop(paste0("Output folder ", dirname(outfile), " does not exist!"))
        }
    }
    if (verbose) {
        cat("Reading in SNP data.\n")
    }
    snps <- read_snp(snpfile)
    snps <- snps %>% dplyr::mutate(relative_pos = 1:dplyr::n())
    if (filter_deam) {
        snps <- snps %>% dplyr::mutate(deam = dplyr::case_when(anc == 
            "C" & der == "T" ~ T, anc == "G" & der == "A" ~ T, 
            T ~ F)) %>% dplyr::filter(!deam) %>% dplyr::select(-deam)
    }
    available_chr <- snps$chr %>% unique()
    if (!is.null(chromosomes)) {
        if (!any(chromosomes %in% available_chr)) {
            stop(paste0("Chromosome names do not match SNP file:\n", 
                "Requested: ", paste0(chromosomes, collapse = ", "), 
                "\n", "Available: ", paste0(available_chr, collapse = ", "), 
                "\n"))
        }
        else {
            snps <- snps %>% dplyr::filter(chr %in% chromosomes)
        }
    }
    if (!is.logical(byPop)) {
        stop("byPop must be logical (TRUE/FALSE)")
    }
    if (verbose) {
        cat(paste0("Analysing chromosomes:\n", paste0(unique(snps$chr), 
            collapse = ", "), "\n"))
    }
    if (is.null(pop_pattern)) {
        ind <- ind_raw %>% dplyr::mutate(row_number = 1:dplyr::n())
    }
    else {
        ind <- ind_raw %>% dplyr::mutate(row_number = 1:dplyr::n()) %>% 
            dplyr::filter(pop %in% pop_pattern)
    }
    n_ind <- nrow(ind)
    if (n_ind <= 1) {
        if (is.null(pop_pattern)) {
            stop(paste0(n_ind, " individuals in final, filtered ind file. Check ind file."))
        }
        else {
            stop(paste0(n_ind, " individuals in final, filtered ind file. Check ind file, and check pop_pattern."))
        }
    }
    geno_list <- vector(mode = "list", length = n_ind)
    if (verbose) {
        cat("Starting to read in genotype data.\n")
    }
    pb1 <- utils::txtProgressBar(min = 1, max = length(geno_list), 
        initial = 1, style = 3)
    dt <- data.table::fread(genofile, col.names = "snp", colClasses = "character")
    for (i in 1:n_ind) {
        dt[, `:=`(digit, stringr::str_sub(snp, ind$row_number[i], 
            ind$row_number[i])), ]
        geno_list[[i]] <- (dt[[2]] %>% dplyr::na_if("9"))[snps$relative_pos]
        utils::setTxtProgressBar(pb1, i)
    }
    if (verbose) {
        cat("  Complete.\n\n")
        cat("Starting to compare genotypes and calculate PMR.\n")
    }
    if (!byPop) {
        pop_group <- rep("Merged", nrow(ind))
    }
    else {
        pop_group <- ind %>% dplyr::group_by(pop) %>% dplyr::mutate(N = dplyr::n()) %>% 
            dplyr::ungroup() %>% dplyr::mutate(pop_group = ifelse(N < 
            minMerge, "Merged", pop)) %>% dplyr::pull(pop_group)
    }
    if (min(table(pop_group)) == 1) {
        stop("Even after possibly merging, smallest population group size is 1.")
    }
    else {
        POPS <- pop_group %>% base::unique() %>% base::sort()
        if (byPop) {
            if (sum(pop_group == "Merged") > 0) {
                ind %>% dplyr::slice(which(pop_group == "Merged")) %>% 
                  dplyr::pull(pop) %>% base::unique() %>% base::sort() %>% 
                  base::paste0(collapse = "/") %>% base::sprintf("Merged group made from: %s\n", 
                  .) %>% message()
            }
        }
    }
    for (k in 1:length(POPS)) {
        if (k == 1) {
            comp_list <- base::vector(mode = "list", length = length(POPS))
        }
        currInds <- ind %>% dplyr::filter(pop_group == POPS[k])
        n_ind <- nrow(currInds)
        if (byPop) {
            sprintf("Comparing within population: %s (%i comparisons)", 
                POPS[k], choose(n_ind, 2)) %>% message()
        }
        counter <- 1
        out_tibble <- tibble::tibble(pair = rep("x", choose(n_ind, 
            2)), nsnps = 0, mismatch = 0, pmr = 0)
        pb2 <- utils::txtProgressBar(min = 0, max = choose(n_ind, 
            2), initial = 1, style = 3)
        for (i in 1:(n_ind - 1)) {
            ind_i <- currInds$ind[i]
            vi <- geno_list[[i]]
            for (j in (i + 1):n_ind) {
                ind_j <- currInds$ind[j]
                vj <- geno_list[[j]]
                snps_filtered_prefiltered <- snps %>% dplyr::mutate(row_num = 1:dplyr::n()) %>% 
                  dplyr::filter(!is.na(vi), !is.na(vj))
                if (nrow(snps_filtered_prefiltered) > 0) {
                  snps_filtered <- snps_filtered_prefiltered
                  out_tibble$pair[counter] <- paste0(ind_i, " - ", 
                    ind_j)
                  out_tibble$nsnps[counter] <- (!is.na(vi[snps_filtered$row_num]) & 
                    !is.na(vj[snps_filtered$row_num])) %>% sum()
                  out_tibble$mismatch[counter] <- (vi[snps_filtered$row_num] != 
                    vj[snps_filtered$row_num]) %>% sum(na.rm = T)
                  out_tibble$pmr[counter] <- out_tibble$mismatch[counter]/out_tibble$nsnps[counter]
                }
                else {
                  out_tibble$pair[counter] <- paste0(ind_i, " - ", 
                    ind_j)
                  out_tibble$nsnps[counter] <- 0
                  out_tibble$mismatch[counter] <- 0
                  out_tibble$pmr[counter] <- NA
                }
                counter <- counter + 1
                utils::setTxtProgressBar(pb2, counter)
            }
        }
        comp_list[[k]] <- out_tibble
    }
    OUT_TIBBLE <- comp_list %>% dplyr::bind_rows()
    if (!is.null(outfile)) {
        if (verbose) {
            cat(sprintf("\n\nSaving processed data to %s.\n", 
                outfile))
        }
        readr::write_delim(OUT_TIBBLE, file = outfile, delim = "\t")
    }
    if (verbose) {
        cat("\nComplete.\n\n")
    }
    OUT_TIBBLE %>% return()
}


# analyze with https://jonotuke.github.io/BREADR/

ind_path <- "analysis/pseudo_haps/highcov.ind"
snp_path <- "analysis/pseudo_haps/highcov.snp"
geno_path <- "analysis/pseudo_haps/highcov.geno"

counts_highcov <- processEigenstratNoFilter(
  indfile = ind_path,
  snpfile = snp_path,
  genofile = geno_path,
  outfile = "analysis/pseudo_haps/highcov.tsv",
  filter_length = 1
)
counts_highcov

relatedness_highcov <- callRelatedness(counts_highcov,
                          class_prior = c(0.05, 0.1, 0.1, 0.75))
relatedness_highcov

plotLOAF(relatedness_highcov, N=300)

plotSLICE(relatedness_highcov, row = 5, labels = c('A', 'B'))


test_degree(relatedness_highcov, row = 5, degree = 3, verbose = TRUE)


BREADR::priorSensitivity(
  in_BREADR = relatedness_highcov,
  row = 10,
  degrees = c(1,2,3,4),
  grid_space = 0.05,
  maxPrior = 0.95
)


##################
# are my related individuals the same as the ones in the paper

fullRel <- read.csv("analysis/pseudo_haps/relatednessNewData.csv", header=T)

fullRel$Indv1 <- gsub("Bricei", "Bric", fullRel$Indv1)
fullRel$Indv2 <- gsub("Bricei", "Bric", fullRel$Indv2)

# use rab

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






#---------------------------------------------------------------------
#---------------------------------------------------------------------
#---------------------------------------------------------------------

# analyze subsampled indivs:

#highcov_Bric001_subsampled_only

#highcov_Bric001_subsampled_merged



ind_path <- "analysis/pseudo_haps/highcov_Bric002_subsampled_merged.ind"
snp_path <- "analysis/pseudo_haps/highcov_Bric002_subsampled_merged.snp"
geno_path <- "analysis/pseudo_haps/highcov_Bric002_subsampled_merged.geno"

counts_highcov <- processEigenstratNoFilter(
  indfile = ind_path,
  snpfile = snp_path,
  genofile = geno_path,
  outfile = "analysis/pseudo_haps/highcov_Bric002_subsampled_merged.tsv"
)
counts_highcov

relatedness_highcov <- callRelatedness(counts_highcov)

relatedness_highcov_orig <- relatedness_highcov

ggplot(relatedness_highcov_orig, aes(x = pmr)) +
  geom_histogram(alpha = 0.7, position = "identity", bins=60) +
  theme_classic(base_size = 14)

relatedness_highcov <- relatedness_highcov |>
  separate(pair, into = c("id1", "id2"), sep = " - ", remove = FALSE) |>
  filter(!(str_detect(id1, "rep") & str_detect(id2, "rep"))) |>
  select(-id1, -id2)

relatedness_highcov <- relatedness_highcov |>
  arrange(pmr)

head(relatedness_highcov$pair)

relatedness_highcov <- relatedness_highcov |>
  separate(pair, into = c("id1", "id2"), sep = " - ", remove = FALSE) |>
  mutate(comparison_type = case_when(
    str_detect(id1, "rep") | str_detect(id2, "rep") &
      (str_detect(id1, "^Bric002$") | str_detect(id2, "^Bric001$")) ~ "subsampled vs highCov Self",
    str_detect(id1, "rep") | str_detect(id2, "rep") ~ "subsampled vs other",
    TRUE ~ "zhighCov only"
  )) |>
  select(-id1, -id2)

table(relatedness_highcov$comparison_type)

ggplot(relatedness_highcov, aes(x = pmr, fill = comparison_type)) +
  geom_histogram(alpha = 0.7, position = "identity", bins=300) +
  theme_classic(base_size = 14)

plotLOAF(relatedness_highcov, N=7688)

plotSLICE(relatedness_highcov, row = 5, labels = c('A', 'B'))


test_degree(relatedness_highcov, row = 5, degree = 3, verbose = TRUE)


BREADR::priorSensitivity(
  in_BREADR = relatedness_highcov,
  row = 10,
  degrees = c(1,2,3,4),
  grid_space = 0.05,
  maxPrior = 0.95
)



# compare to python pmr values:

# Python output
pmr_python <- read.csv("analysis/pseudo_haps/pmr_merged_bric002_python.csv")

# assumes df_pmr_merged had columns like iid1, iid2, mms, snps, pmr —
# adjust names below if yours differ
pmr_python <- pmr_python |>
  mutate(id1 = pmin(iid1, iid2),
         id2 = pmax(iid1, iid2)) |>
  select(id1, id2, nsnps_py = snps, mismatch_py = mms, pmr_py = pmr)

counts_highcov_nofilter <- processEigenstratNoFilter(
  indfile = ind_path,
  snpfile = snp_path,
  genofile = geno_path,
  outfile = NULL
)

pmr_r <- counts_highcov_nofilter |>
  separate(pair, into = c("Indv1", "Indv2"), sep = " - ", remove = FALSE) |>
  mutate(id1 = pmin(Indv1, Indv2),
         id2 = pmax(Indv1, Indv2)) |>
  select(id1, id2, nsnps_r = nsnps, mismatch_r = mismatch, pmr_r = pmr)

# merge and compare
comparison <- pmr_r |>
  inner_join(pmr_python, by = c("id1", "id2")) |>
  mutate(
    pmr_diff = pmr_r - pmr_py,
    nsnps_diff = nsnps_r - nsnps_py
  )

# sanity checks
summary(comparison$pmr_diff)
summary(comparison$nsnps_diff)

# flag any pairs that disagree substantially
comparison |>
  filter(abs(pmr_diff) > 0.01) |>
  arrange(desc(abs(pmr_diff)))

# specifically pull out the self-vs-replicate pairs you were debugging
comparison |>
  filter(str_detect(id1, "rep") | str_detect(id2, "rep")) |>
  arrange(pmr_r)


threshold <- 2

comparison_filtered <- comparison |>
  filter(nsnps_r > threshold, nsnps_py > threshold)|>
  arrange(desc(abs(pmr_diff)))
comparison_filtered
comparison_filtered |>
  filter(abs(pmr_diff) > 0.01) |>
  arrange(desc(abs(pmr_diff)))

hist(comparison$pmr_diff, xlim=c(-0.1, 0.1))







####### ------------------------
ind_path <- "analysis/pseudo_haps/all_samples_merged.ind"
snp_path <- "analysis/pseudo_haps/all_samples_merged.snp"
geno_path <- "analysis/pseudo_haps/all_samples_merged.geno"

counts_highcov <- processEigenstratNoFilter(
  indfile = ind_path,
  snpfile = snp_path,
  genofile = geno_path,
  outfile = "analysis/pseudo_haps/all_merged.tsv",
  filter_length = 1
)
counts_highcov

relatedness_all <- callRelatedness(counts_highcov)



# lists of who's who, straight from your Python sample lists
highcov_ids <- c("Bric001","Bric002","Bric003","Bric004","Bric005","Bric007","Bric008",
                  "Bric011","Bric015","Bric016","Bric021","Bric022","Bric023","Bric024",
                  "Bric031","Bric033","Bric034","Bric035","Bric037","Bric039","Bric040",
                  "Bric042","Bric043","Bric044","Bric006","Bric026","Bric032")
lowcov_ids <- c("Bric045","Bric046","Bric-eDNA03","Bric-eDNA04","Bric-eDNA05","Bric-eDNA06",
                "Bric-eDNA07","Bric-eDNA08","Bric-eDNA09","Bric-eDNA10","Bric-eDNA11",
                "Bric-eDNA12","Bric-eDNA13","Bric-tag01","Bric-tag02","Bric-tag03")

relatedness_all <- relatedness_all |>
  separate(pair, into = c("id1", "id2"), sep = " - ", remove = FALSE) |>
  mutate(comparison_type = case_when(
    id1 %in% lowcov_ids & id2 %in% lowcov_ids ~ "lowCov vs lowCov",
    id1 %in% lowcov_ids | id2 %in% lowcov_ids ~ "lowCov vs highCov",
    TRUE ~ "highCov only"
  )) |>
  select(-id1, -id2)

table(relatedness_all$comparison_type)

ggplot(relatedness_all, aes(x = pmr, fill = comparison_type)) +
  geom_histogram(alpha = 0.7, position = "identity", bins = 300) +
  theme_classic(base_size = 14)


relatedness_all_filtered <- relatedness_all |>
  filter(comparison_type != "lowCov vs lowCov")|>
  filter(comparison_type != "highCov only")

table(relatedness_all_filtered$comparison_type)

ggplot(relatedness_all_filtered, aes(x = pmr, fill = comparison_type)) +
  geom_histogram(alpha = 0.7, bins = 300) +
  theme_classic(base_size = 14)


head(relatedness_all_filtered)

out <- relatedness_all_filtered |>
  arrange(pmr)

out[1:20,1:7]

# everything is a twin with 006, which is the reference genome. Reference bias?


# drop 006

relatedness_all_filtered2 <- relatedness_all_filtered |>
  filter(!str_detect(pair, "Bric006"))

ggplot(relatedness_all_filtered2, aes(x = pmr, fill = comparison_type)) +
  geom_histogram(alpha = 0.7, bins = 300) +
  theme_classic(base_size = 14)

out <- relatedness_all_filtered2 |>
  arrange(pmr)

print(out[1:30,1:7], n=30)


ggplot(relatedness_all_filtered3, aes(x = pmr, fill = comparison_type)) +
  geom_histogram(alpha = 0.7, bins = 300) +
  theme_classic(base_size = 14)

out <- relatedness_all_filtered3 |>
  arrange(pmr)

out[1:20,1:7]
