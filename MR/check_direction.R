#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
output <- args[1]
input1 <- args[2]
input2 <- args[3]

# print args
message("output = ", output, "\ninput1 = ", input1, "\ninput2 = ", input2)

mr_df <- fread(input1) %>%
    filter(is.na(pval) == FALSE)
#mr_df %>% arrange(pval) %>% select(outcome, exposure, pval) %>% head(20)

# check the inverse 
other_df <- fread(input2) %>%
    filter(is.na(pval) == FALSE)

# make new df
checked_df <- data.frame()
for (i in seq_len(nrow(mr_df))) {
    row <- mr_df[i, ]
    other_row <- other_df %>% filter(outcome == row$exposure, exposure == row$outcome)
    if (nrow(other_row) != 0) {
        if (other_row$pval < 0.05) {  # always gonna be the case if it exists in the file
            message("\nFound inverse: ", row$outcome, " on ", row$exposure, 
                "\n\tThis pval: \t", row$pval, 
                "\n\tOther pval \t", other_row$pval)
            next
        }
    }
    checked_df <- rbind(checked_df, row)
}

#out <- paste0(outdir, "/", pheno1, "_on_", pheno2, "/all_mr_checked.tsv")
fwrite(checked_df, output, sep = "\t")
message(output)
