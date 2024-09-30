#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

root <- "/vol/projects/CIIM/meta_cQTL/out"
cohort <- "2000HIV-EU-discovery"
cov <- "main"
pheno_e <- "proteins"    # exposure
pheno_o <- "cytokines"   # outcome
input <- paste0(root, "/", cohort, "/mr/", cov, "/", pheno_e, "_on_", pheno_o, "/all_mr_checked.tsv")
outdir <- paste0(root, "/", cohort, "/mr/", cov, "/", pheno_e, "_on_", pheno_o, "/per_outcome")

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
outdir <- args[2]
pheno_o <- args[3]

# make outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message(input)
infile <- fread(input)

# load outcome file to know all the possible outcomes
outcomes <- fread(paste0(root, "/", cohort, "/", pheno_o, "/phenotype/normalised.tsv")) %>%
    pull(V1) %>% 
    unique() %>% 
    sort()

# split infile per outcome and save to file
o <- outcomes[1]
for (o in outcomes) {
    outfile <- paste0(outdir, "/", o, ".tsv")
    message(outfile)
    infile %>% filter(outcome == o) %>% fwrite(outfile, sep = "\t")
}

