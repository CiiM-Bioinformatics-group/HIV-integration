#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
input <- as.list(args)[-c(1:2)] # everything except 1st arg

message("outfile: ", outfile, 
        "\ninput: ", paste(input, collapse = "\n"))

# create outdir
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

# read qtl results
qtl <- fread(cmd = paste0("awk 'NR == 1 || FNR > 1' ", paste(input, collapse = " ")))
names(qtl) <- c("SNP", "gene", "beta", "t-stat", "p")

# calculate lambda
lambda <- qtl %>% summarise("lam" = median(qchisq(1 - na.omit(qtl$p), 1)) / qchisq(0.5, 1))

# save lambda
fwrite(lambda, outfile, col.names = FALSE)
message("done", outfile)