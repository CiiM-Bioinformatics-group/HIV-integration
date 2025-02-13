#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

message("infile ", infile,
        "\noutfile ", outfile
)

# make plot dir
dir.create(dirname(outfile), showWarnings = FALSE)

df <- fread(infile, header = TRUE) %>% column_to_rownames(names(.)[1])
nrow(df)

message("INT normalization")
int_norm <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
df <- apply(df, 1, int_norm) %>% t() %>% as.data.frame()


fwrite(df, outfile, row.names = TRUE, sep = "\t")
message(outfile)
