#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

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
