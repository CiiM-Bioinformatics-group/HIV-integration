#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

message("infile: ", infile,
        "\noutfile: ", outfile)

# from javi
meff_li <- function(mat) {
    #mat: a matrix having the variables (traits, cytokines, metabolites...) as columns
    #Formula 5 taken from J Li and L Ji Heredity 2005

    #Calculate correlation matrix
    # if mat contains any NA
    if (mat %>% is.na %>% any) {
        corr_matrix <- cor(mat, method = "spearman", use = "pairwise.complete")
    } else {
        corr_matrix <- cor(mat, method = "spearman")
    }
    #Any correlation of NA (no matching values) is set to 0
    corr_matrix[is.na(corr_matrix)] <- 0
    #Eigenvalues for correlation matrix
    evals <- eigen(t(corr_matrix), symmetric = TRUE)$values
    #Formula from the paper to calculate the number of effective tests (meff)
    f <- function(x) {
        ifelse(x < 0, 0, ifelse(x >= 1, 1, 0)) + (x - floor(x))
    }
    meff <- sum(f(evals))
    #Alpha value calculation
    alpha <- 1 - (1 - 0.05) ^(1 / (meff))
    alpha_gw <- alpha / 1e6

    return(list("meff" = meff, "alpha" = alpha, "alpha_gw" = alpha_gw))
}

# load in phenotype
pheno_infile <- fread(infile) %>% column_to_rownames("V1")
pheno_t <- pheno_infile %>% t() %>% as.data.frame()

# calc
ret <- meff_li(pheno_t)

# write results
fwrite(ret, outfile, sep = "\t")
message(outfile)
