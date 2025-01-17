#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

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
