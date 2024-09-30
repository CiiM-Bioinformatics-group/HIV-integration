#!/usr/bin/env Rscript

# this script takes 1-e5 snps and performs very strong clumping to get independent snps

library(tidyverse)
library(data.table)

library(future)
library(future.apply)

library(ieugwasr)
library(plinkbinr)
#library(LDlinkR)

options(future.globals.maxSize = 8000 * 1024^2)

root <- "/vol/projects/CIIM/meta_cQTL/out"
cohort <- "2000HIV-EU-discovery"
cov <- "main"
pheno <- "expression"  
bfile <- paste0(root, "/", cohort, "/genotype/allchr")
pval <- "1e-5" 
outpath <- paste0(root, "/", cohort, "/mr/", cov, "/clump/", pheno)
val_cohort <- "2000HIV-EU-validation"
val_pval_threshold <- 0.05
chr_num <- 21

args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
cohort <- args[2]
cov <- args[3]
pheno <- args[4]
bfile <- args[5]
pval <- args[6]
outpath <- args[7]
val_cohort <- args[8]
val_pval_threshold <- as.numeric(args[9])
chr_num <- as.numeric(args[10])

dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# gather the effect allele frequency for each snp ;; using future
tmp_path <- file.path(root, cohort, "imputation", "freq", paste0("chr", chr_num, ".afreq"))
all_eaf_dat <- fread(tmp_path, header = TRUE) %>%
        separate(ID, c("other", "snp"), sep = ";", fill = "right") %>% 
        mutate(MAF = ifelse(ALT_FREQS < 0.5, ALT_FREQS, 1 - ALT_FREQS)) %>%
        #filter(MAF >= 0.05) %>%
        dplyr::rename(eaf = ALT_FREQS) %>%
        dplyr::select(snp, eaf)
message("Processed chr ", chr_num)

# check memory usage
print(object.size(x = lapply(ls(), get)), units = "Mb")

cmd <- paste0("awk -F : '{if ($1 == \"chr", chr_num, "\") print}' ", root, "/", cohort, "/", pheno, "/mapping/main_", pval, ".tsv")
exposure_snps <- fread(cmd = cmd, header = TRUE) %>%
    setNames(c("SNP", "gene", "beta", "t-stat", "p-value")) %>%
    separate(SNP, c("chr", "location", "ref", "alt", "snp"), sep = "[:;]", fill = "right", remove = FALSE) %>%
    dplyr::filter(chr == chr_num) %>%
    dplyr::rename("pval" = `p-value`) %>%
    dplyr::filter(!is.na(snp)) %>%  #remove any where snp = NA
    mutate(se = abs(beta) / `t-stat`, location = as.numeric(location)) %>%
    mutate(pheno = pheno)
nrow(exposure_snps) # 

# check memory usage
print(object.size(x = lapply(ls(), get)), units = "Mb")

exposure_snps <- exposure_snps %>% inner_join(all_eaf_dat, by = "snp")
nrow(exposure_snps) # 433608 -- why?? -- oh right im filtering MAF >= 0.05

# check memory usage
print(object.size(x = lapply(ls(), get)), units = "Mb")

#unload all_eaf_dat
rm(all_eaf_dat)

# check memory usage
print(object.size(x = lapply(ls(), get)), units = "Mb")

# remove snps not nominal significant in validation
# get all exposure snps for this chromosome
validation_dat <- fread(paste0(root, "/", val_cohort, "/", pheno, "/mapping/", cov, "/nominal/chr", chr_num, ".tsv"))
nrow(validation_dat) # 105825306

if (nrow(validation_dat) == 0) {
    return(data.frame())
}
validation_dat <- validation_dat %>%
    dplyr::select(SNP, gene, "val_pval" = `p-value`)

validated_snps <- exposure_snps %>%
    filter(chr == paste0("chr", chr_num)) %>%
    inner_join(validation_dat, by = c("SNP", "gene"))

# check memory usage
print(object.size(x = lapply(ls(), get)), units = "Mb")

# get all unique ones from phenotype file to not miss any 
traits <- fread(paste0(root, "/", cohort, "/", pheno, "/phenotype/normalised.tsv"), header = TRUE) %>%
    pull(V1) %>%
    unique() %>%
    sort()

# split on gene (trait) and clump
plan(multisession)
trait1 <- "CRYZL1_Cardiometabolic_II"
trait1 <- exposure_snps$gene[1]
future_lapply(traits, future.seed = TRUE, function(trait1) {
    # get snps for each trait
    exposure_dat_raw <- validated_snps %>%
        filter(gene == trait1) %>% 
        unique() #sometimes theres duplicate rows, idk why yet, eg. expression ENSG00000249601.2
    nrow(exposure_dat_raw)

    if (nrow(exposure_dat_raw > 0)) {
        # clump
        exposure_dat_clump <- exposure_dat_raw %>% 
            dplyr::rename(rsid = "SNP") %>%
            ld_clump(clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, 
                bfile = bfile, plink_bin = plinkbinr::get_plink_exe()) #genetics.binaRies::get_plink_binary())
        exposure_dat <- exposure_dat_clump %>%
            dplyr::rename(ID = rsid, SNP = snp, pval.exposure = "pval", beta.exposure = "beta",
            chr.exposure = "chr", pos.exposure = "location", effect_allele.exposure = "ref",
            other_allele.exposure = "alt", se.exposure = "se", eaf.exposure = "eaf", 
            exposure = "gene", id.exposure = "pheno") %>%
            arrange(ID)
        #if (nrow(exposure_dat) >= 3) {
        # write to file
        fwrite(exposure_dat, file.path(outpath, paste0(trait1, ".tsv")), sep = "\t")
        #}
    } else {
        # make empty file
        fwrite(data.frame(), file.path(outpath, paste0(trait1, ".tsv")), sep = "\t")
    }
    trait1
})
