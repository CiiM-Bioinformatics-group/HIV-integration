#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

outdir <- "/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/mr/main"
pheno1 <- "expression"    # exposure
pheno2 <- "cytokines"   # outcome

output <- paste0(outdir, "/", pheno1, "_on_", pheno2, "/all_mr_checked_pval.tsv")
input <- paste0(outdir, "/", pheno1, "_on_", pheno2, "/all_mr_checked.tsv")
pval_file <- paste0(outdir, "/", pheno1, "_on_", pheno2, "/pvalue.txt")

args <- commandArgs(trailingOnly = TRUE)
output <- args[1]
input <- args[2]
pval_file <- args[3]

# print args
message("output = ", output, "\ninput = ", input, "\npval_file = ", pval_file)

df <- fread(input)
nrow(df)

pval_threshold <- fread(pval_file, header = FALSE)$V1

df_filter <- df %>% filter(pval < pval_threshold)
nrow(df_filter)

df_filter %>% fwrite(output, sep = "\t", na = "NA", quote = FALSE, row.names = FALSE)
message(output)
