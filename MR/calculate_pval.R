#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

outdir <- "/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery"
pheno1 <- "expression"    # exposure
pheno2 <- "cytokines"   # outcome
cov <- "main"

input1 <- paste0(outdir, "/", pheno1, "/phenotype/multiple_testing.tsv")
input2 <- paste0(outdir, "/", pheno2, "/phenotype/multiple_testing.tsv")
output <- paste0(outdir, "/mr/", cov, "/", pheno1, "_on_", pheno2, "/pvalue.txt")

args <- commandArgs(trailingOnly = TRUE)
input1 <- args[1]
input2 <- args[2]
output <- args[3]

# print args
message("input1 = ", input1, "\ninput2 = ", input2, "\noutput = ", output)

n_1 <- fread(input1, header = TRUE)$meff
n_1
n_2 <- fread(input2, header = TRUE)$meff
n_2

pval_threshold <- 0.05 / (n_1 * n_2)
pval_threshold
fwrite(pval_threshold %>% as.data.frame(), output, row.names = FALSE, col.names = FALSE)
message(output)
