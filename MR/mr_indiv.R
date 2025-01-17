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

args <- commandArgs(trailingOnly = TRUE)
file1 <- args[1]
file2 <- args[2]
pheno1 <- args[3]
pheno2 <- args[4]
outpath <- args[5]

message("file1 = ", file1, 
        "\nfile2 = ", file2, 
        "\npheno1 = ", pheno1, 
        "\npheno2 = ", pheno2, 
        "\npval_val = ", pval_val, 
        "\noutpath = ", outpath,
        "\ntesting = ", testing,
        "\nn_samples = ", n_samples,
        "\n")

dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
message(outpath)

# clumping is done ahead with genomewide_clumping.R
trait1 <- file1 %>% basename() %>% str_remove("\\.tsv")
message(trait1)

# distinct cause for some eqtl it gives duplicate rows, same snp, beta, everything
exposure_df <- fread(file1) %>% distinct()
# invert eaf.exposure
exposure_df <- exposure_df %>% mutate(eaf.exposure = 1 - eaf.exposure)

if (nrow(exposure_df) == 0) {
    file.create(paste0(outpath, "/all_mr.tsv"))
} else {
    # sometimes lines occur twice but with different pval so distinct() doesn't catch
    exposure_dup <- duplicated(exposure_df %>% dplyr::select(-pval.exposure))
    exposure_raw <- exposure_df[!exposure_dup]
    
    # check pleiotropy-- to be moved to after clumping but now temporarily here
    # load in the whole -5 qtl file

    if (FALSE) {# moved to clumping step
        qtl_file <- paste0(root, "/", cohort, "/", pheno1, "/mapping/main_1e-5.tsv")
        qtls <- fread(qtl_file)
        
        # for snp in exposure
        rets <- data.frame()
        i <- 1
        for (i in 1:nrow(exposure_raw)) {
            message(i)
            row <- exposure_raw[i,]
            # find if SNP has more significant associations
            ret <- qtls %>% dplyr::filter(SNP == row$ID)
            min <- ret %>% dplyr::filter(`p-value` == min(ret$`p-value`))
            min$tot <- nrow(ret)
            rets <- rbind(rets, min)
        }
        print(rets)
        rets <- rets %>% dplyr::select(SNP, tot)
        exposure_dat <- exposure_raw %>% inner_join(rets, by = join_by("ID" == "SNP")) %>% dplyr::filter(tot <= 4)
        message("SNP difference after pleiotropy check: ", nrow(exposure_raw), " -> ", nrow(exposure_dat))
    } else {
        exposure_dat <- exposure_raw 
    }

    outcome_all_raw <- fread(file2) %>% distinct()
    nrow(outcome_all_raw) # 10267
    outcome_all_raw <- outcome_all_raw[!duplicated(outcome_all_raw %>% dplyr::select(-pval))]
    nrow(outcome_all_raw) # 10258

    message("Checking if there are any trait2 with >= 3 snps ...")
    keep2 <- outcome_all_raw %>% group_by(gene) %>% summarise(n = n()) %>% filter(n >= 3) %>% pull(gene)
    length(keep2) # 326 proteins
    # get ones with hivenv in the name
    #keep2[grepl("hivenv", keep2)]

    outcome_all <- outcome_all_raw %>% filter(gene %in% keep2)

    if (nrow(outcome_all) == 0) {
        file.create(paste0(outpath, "/all_mr.tsv"))
    } else {
        dir.create(outpath, showWarnings = FALSE, recursive = TRUE)

        # loop through each unique "gene" in outcome and perform mr per exposure&outcome
        full_sub <- data.frame()
        #failed <- data.frame() %>% mutate(trait1 = character(), trait2 = character(), reason = character())    
        message("Looping through traits from other outcome...")
        if (testing) trait2 <- outcome_all$gene[[1]]  # testing
        if (testing) trait2 <- "pbmc_24h_il1b_hivenv"
        if (testing) trait2 <- "ITGA6_Inflammation"
        if (testing) trait2 <- "ENSG00000203761"
        for (trait2 in unique(outcome_all$gene)) {
            outcome <- outcome_all %>% filter(gene == trait2)
            nrow(outcome) # 1

            outcome_dat <- outcome %>%
                dplyr::rename(pval.outcome = "pval", beta.outcome = "beta",
                chr.outcome = "chr", pos.outcome = "location", effect_allele.outcome = "effect_allele",
                other_allele.outcome = "other_allele", se.outcome = "se", eaf.outcome = "eaf.exposure", 
                outcome = "gene", id.outcome = "pheno")
                
            # sometimes this version works instead:
            #outcome_dat <- outcome %>% rename(pval = "pval.outcome", beta = "beta.outcome", chr = "chr.outcome", 
            #    location = "pos.outcome", effect_allele = "effect_allele.outcome", other_allele = "other_allele.outcome",
            #    se = "se.outcome", eaf.exposure = "eaf.outcome", gene = "outcome", pheno = "id.outcome")

            dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
            dat$samplesize.outcome <- n_samples
            dat$samplesize.exposure <- n_samples

            if (nrow(dat) == 0) {
                #failed <- rbind(failed, data.frame(trait1 = trait1, trait2 = trait2, 
                #    reason = "No snps after harmonising"))
                message("No snps after harmonising, skipping ...")
                next
            }
            # filter out mr_keep = FALSE
            dat <- dat %>% filter(mr_keep == TRUE)
            if (nrow(dat) < 3) {
                #failed <- rbind(failed, data.frame(trait1 = trait1, trait2 = trait2, 
                #    reason = "Less than 3 snps after harmonising"))
                message("Less than 3 snps after harmonising, skipping ...")
                next
            }

            outname <- paste0(outpath, "/", trait2) #ivw_pval %>% format(scientific = TRUE, digits = 3), "_",
            dir.create(outname, showWarnings = FALSE, recursive = TRUE)

            # save data
            fwrite(dat, paste0(outname, "/data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

            # perform MR
            all_res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_simple_median", "mr_weighted_median"))
            res <- all_res %>% filter(method == "Inverse variance weighted") %>%
                mutate(pleio = NA, het = NA, lou = NA)
            ivw_pval <- res %>% pull(pval)

            # if significant, perform sensitivity checks and plots
            if (ivw_pval < 0.05) {
                # single snp  & forest
                s_res <- mr_singlesnp(dat)
                tryCatch({
                    frplot <- mr_forest_plot(s_res)
                    ggsave(plot = frplot[[1]], filename = paste0(outname, "/forest.pdf"), width = 4, height = 5)
                }, error = function(e) {
                    # e.g. Error in `levels<-`(`*tmp*`, value = as.character(levels)): factor level [5] is duplicated
                    # happens when directionality also cannot be calculated
                    # and there are warnings on perfect fits... 
                    # e.g. In summary.lm(mod) : essentially perfect fit: summary may be unreliable
                    message("Error in forest plot ...\n", e)
                })
                
                # Plot MR
                out_plot <- paste0(outname, "/scatter.pdf")
                mrplot <- mr_scatter_plot(all_res, dat)[[1]] + 
                    theme_bw() + 
                    labs(title = paste0("MR ", trait1, " vs ", trait2), 
                        subtitle = paste0("IVW p-value: ", ivw_pval %>% format(scientific = TRUE, digits = 3))
                                        #"\nDirectionality p-value: ", dir$steiger_pval %>% format(scientific = TRUE, digits = 3),
                                        #"\nHeterogeneity p-value: ", het$Q_pval %>% format(scientific = TRUE, digits = 3), 
                                        #"\nPleiotropy p-value: ", pleio$pval %>% format(scientific = TRUE, digits = 3),
                                        #"\nLeave-one-out max p-value: ", max(lou$p) %>% format(scientific = TRUE, digits = 3)
                                        ) +
                        geom_vline(aes(xintercept = 0), linetype = "dashed") + 
                        geom_hline(aes(yintercept = 0), linetype = "dashed")
                ggsave(plot = mrplot, filename = out_plot, width = 10, height = 5)
                message("Saved plot to ", out_plot)

                # sensitivity checks
                pleio <- mr_pleiotropy_test(dat)
                res$pleio <- pleio$pval
                het <- mr_heterogeneity(dat, method_list = "mr_ivw")
                res$het <- het$Q_pval
                lou <- mr_leaveoneout(dat)
                if (!is.null(lou)) {
                    if (any(is.na(lou$p))) {
                        res$lou <- NA
                    } else {
                        res$lou <- max(lou$p)
                    }
                }
            }

            full_sub <- rbind(full_sub, res)
        }

        outfile <- paste0(outpath, "/all_mr.tsv")
        fwrite(full_sub, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
        message("Saved results to ", outfile)
        #full <- rbind(full, full_sub)

        #if (nrow(failed) > 0) {
        #    outfile <- paste0(outpath, "/all_failed.tsv")
        #    fwrite(failed, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
        #    message("Saved failed to ", outfile)
        #}
        #if (length(dir(outpath)) == 0) {
        #    unlink(outpath, recursive = TRUE)
        #}
    }
}
message("done ")
