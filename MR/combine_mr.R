#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
pheno1 <- args[3]
pheno2 <- args[4]

message(indir)

all_mr <- data.frame()
for (subdir in list.dirs(indir, recursive = FALSE)) {
    files <- list.files(subdir, recursive = FALSE)

    #} else if (length(files) > 1 && "all_mr.tsv" %in% files) { 
    # combine all all_mr.tsv into one
    message("Combining ", subdir)
    tmp_mr <- files %>% 
        str_subset("all_mr.tsv") %>% 
        map_df(~fread(paste0(subdir, "/", .x)))
    all_mr <- rbind(all_mr, tmp_mr)
    #}
}

# Benjamini-Hochberg (FDR) correction
all_mr$pval_adj <-  p.adjust(all_mr$pval, method = "BH")

if (nrow(all_mr) != 0) {
    all_mr <- all_mr %>% arrange(pval)
} else {
    message("No mr found")
}

sign_mr <- all_mr %>% filter(pval < 0.05)
dim(sign_mr) # 1695   13

sign_mr_check <- sign_mr %>% filter(het > 0.05 & pleio > 0.05 & lou < 0.05)
dim(sign_mr_check) # 135   13

sign_mr_adj <- sign_mr_check %>% filter(pval_adj < 0.05)
dim(sign_mr_adj) # 0 13

fwrite(all_mr, file.path(outdir, "all_mr.tsv"), sep = "\t")
message("Done ", file.path(outdir, "all_mr.tsv"))

fwrite(sign_mr_check, file.path(outdir, "significant_mr.tsv"), sep = "\t")
message("Done ", file.path(outdir, "significant_mr.tsv"))

fwrite(sign_mr_adj, file.path(outdir, "significant_mr_adj.tsv"), sep = "\t")
message("Done ", file.path(outdir, "significant_mr_adj.tsv"))


