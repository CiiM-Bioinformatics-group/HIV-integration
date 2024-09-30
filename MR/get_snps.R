#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(ggplot2)

library(future)
library(future.apply)

library(TwoSampleMR) # mr
library(ieugwasr)
library(plinkbinr)
#library(LDlinkR)

options(future.globals.maxSize = 8000 * 1024^2)

root <- "/vol/projects/CIIM/meta_cQTL/out"
cohort <- "2000HIV-EU-discovery"
cov <- "main"
pheno1 <- "cytokines"    # exposure
pheno2 <- "expression"   # outcome
#file <- paste0(root, "/", cohort, "/mr/", cov, "/clump/", pheno1, "/A1BG_Inflammation_II.tsv")
#file <- paste0(root, "/", cohort, "/mr/", cov, "/clump/", pheno1, "/Panel3_3116_Naive_B_cells_CD307d_.tsv")
file <- paste0(root, "/", cohort, "/mr/", cov, "/clump/", pheno1, "/pbmc_7d_il5_spneu.tsv")
outpath <- paste0(root, "/", cohort, "/mr/", cov, "/", pheno1, "_on_", pheno2)

file = "/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/mr/main/clump/cytokines/pbmc_24h_mcp1_hivenv.tsv"
pheno1 = "cytokines"
pheno2 = "metabolites"
outpath = "/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/mr/main/cytokines_on_metabolites/data"

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
pheno1 <- args[2]
pheno2 <- args[3]
outpath <- args[4]
root <- args[5]
cohort <- args[6]
cov <- args[7]

message("file = ", file, "\npheno1 = ", pheno1, "\npheno2 = ", pheno2, "\noutpath = ", outpath)

dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
message(outpath)

# clumping is done ahead with genomewide_clumping.R
trait1 <- file %>% basename() %>% str_remove("\\.tsv")
message(trait1)

exposure_dat <- fread(file)

if (nrow(exposure_dat) == 0) {
    message("No SNPs for ", trait1)
    outcome_all <- data.frame()
} else {
    # Get unique chr.exposure values from exposure_dat
    unique_chr_exposure <- unique(exposure_dat$chr.exposure)

    # Function to process data for each unique chr.exposure
    process_data <- function(chr_exposure) {
        if (pheno2 == "cell_count") {
            cmd <- paste0("awk -F '\t' 'NR == 1 || FNR > 1 && ", 
                    paste0("$1 == \"", exposure_dat$ID[exposure_dat$chr.exposure == chr_exposure], collapse = "\" || "),
                    "\"' ", root, "/", cohort, "/", pheno2, "/mapping-chr/", cov, "/", chr_exposure, "/*.tsv")
        } else if (pheno2 == "expression") {
            cmd <- paste0("awk -F '\t' 'NR == 1 || FNR > 1 && ", 
                    paste0("$1 == \"", exposure_dat$ID[exposure_dat$chr.exposure == chr_exposure], collapse = "\" || "),
                    "\"' ", root, "/", cohort, "/", pheno2, "/mapping/", cov, "/", chr_exposure, "/*.tsv")
        } else {
            cmd <- paste0("awk -F '\t' 'NR == 1 || ", 
                    paste0("$1 == \"", exposure_dat$ID[exposure_dat$chr.exposure == chr_exposure], collapse = "\" || "),
                    "\"' ", root, "/", cohort, "/", pheno2, "/mapping/", cov, "/", chr_exposure, ".tsv")
        }
        tmp <- fread(cmd = cmd)
        message("Done with ", chr_exposure)
        tmp
    }

    # then get all the same snps from the other data 
    message("Getting all snps from other data ...")
    plan(multisession)
    outcome_all_list <- future_lapply(unique_chr_exposure, process_data)
    outcome_all_raw <- do.call(rbind, outcome_all_list)
    nrow(outcome_all_raw) # 661

    # fix some things
    message("Fixing some things ...")
    outcome_all <- outcome_all_raw %>%
        separate(SNP, c("chr", "location", "effect_allele", "other_allele", "SNP"), sep = "[:;]") %>%
        filter(!is.na(SNP)) %>%  #remove any where snp = NA
        rename(pval = "p-value") %>%
        mutate(se = abs(beta) / `t-stat`, location = as.numeric(location)) %>%
        inner_join(exposure_dat %>% 
                    select(SNP, eaf.exposure), by = "SNP") %>%
        mutate(pheno = pheno2)
}
fwrite(outcome_all, paste0(outpath, "/", trait1, "_snps.tsv"), sep = "\t")
