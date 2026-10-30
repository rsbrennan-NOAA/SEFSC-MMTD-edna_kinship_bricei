---
title: "Rice's whale eDNA WGS"
subtitle: "Testing Whole genome data from environmental and non-biopsy samples"
author: "Reid Brennan"
date: today
date-format: "YYYY-MM-DD"
format:
  html:
    theme: cosmo
    toc: true
    toc-depth: 3
    toc-location: left
    embed-resources: true
    code-fold: true
    code-summary: "Show code"
    code-tools: true
    df-print: paged
    fig-width: 8
    fig-height: 5
    lightbox: true
    anchor-sections: true
    smooth-scroll: true
execute:
  echo: true
  warning: false
  message: false
  cache: true
---


# Rice's whale eDNA wgs

Goal: can we use non-biopsy samples to obtain whole genome data from Rice's whales? 
Samples being tested:
- Tube juice
- fluke prints
- Bones

The hope is that we can obtain 

1. Mitochondrial genomes for every samples
2. Low coverage WGS. At a minimum identify duplicate individuals. Ideally with good enough data to use for population genetics. 


## Samples and methods

## Initial results: Mapping

### Number of reads and mapping success

```{r}
#| label: load-data

library(tidyverse)

qc <- read_tsv("qc_report.tsv") |>
  rename(genome_cov = `% genome_cov`) |>
  mutate(
    pct_genome_cov = genome_cov * 100,
    sample_type = case_when(
      sample %in% c("Bric026", "Bric032") ~ "Biopsy",
      sample %in% c("Bric045", "Bric046") ~ "Bone",
      str_detect(sample, "^Bric_tag") ~ "Tube juice",
      str_detect(sample, "^Bric_eDNA") ~ "Fluke print"
    ),
    sample_type = factor(
      sample_type,
      levels = c("Biopsy", "Bone", "Tube juice", "Fluke print")
    )
  )
```

### Reference Bias

### what are the falsely aligned reads?

## Successful Samples

### Tube Juice and Bones



