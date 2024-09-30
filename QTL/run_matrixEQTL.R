#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(bnstruct)
  library(patchwork)
  library(MatrixEQTL)
  library(parallel)
})

base <- "/vol/projects/CIIM/meta_cQTL"
cohort <- "2000HIV-EU-discovery"
cohort_out <- "2000HIV-EU-discovery_old"
pheno_type <- "cytokines"
cov_interact <- "main"
p_val <- 1
cov_main_unsplit <- "age,sex,bmi,lockdown,hiv_duration,season_sin,season_cos"


# proteins or metabolites
pheno_type <- "proteins"
out_mapping <- paste(base, "out", cohort_out, pheno_type, "mapping", cov_interact, "chr22", sep = "/")

# cytokines
pheno_type <- "cytokines"
out_mapping <- paste(base, "out", cohort_out, pheno_type, "mapping-chr", cov_interact, "chr22", sep = "/")


# other paths
snp_path <- paste(base, "out", cohort, "genotype/dosage/chr3.txt", sep = "/")
cov_path <- paste(base, "data", cohort, pheno_type, "covariates.tsv", sep = "/")
pheno_path <- paste(base, "out", cohort, pheno_type, "phenotype/normalised.tsv", sep = "/")


# Overrides the manual paths above which are just used for testing
args <- commandArgs(trailingOnly = TRUE)
snp_path <- args[1]
cov_path <- args[2]
pheno_path <- args[3]
out_mapping <- args[4]
p_val <- as.numeric(args[5])
cov_main_unsplit <- args[6]
cov_interact <- args[7]
pheno_type <- args[8]

# print for checking
message("snp_path ", snp_path,
        "\ncov_path ", cov_path,
        "\npheno_path ", pheno_path,
        "\nout_mapping ", out_mapping,
        "\np_val ", p_val,
        "\ncov_main ", cov_main_unsplit,
        "\ncov_interact ", cov_interact,
        "\npheno_type ", pheno_type)

# functions
create_sliced_data <- function(mat) {
  # create data type required for matrixeqtl
  slice <- SlicedData$new()
  slice$fileOmitCharacters <- "NA"
  slice$fileDelimiter <- "\t"
  slice$fileSkipRows <- 0
  slice$fileSkipColumns <- 0
  slice$fileSliceSize <- 200
  slice$CreateFromMatrix(as.matrix(mat))
  return(slice)
}

perform_mapping <- function(model, pheno_name, this_pheno, this_lb, this_cov, outfile) {
  snps <- create_sliced_data(this_lb)
  pheno <- create_sliced_data(this_pheno)
  cvrt <- create_sliced_data(this_cov)
  quest_qtl <- Matrix_eQTL_engine(
    snps = snps,
    gene = pheno,
    cvrt = cvrt,
    output_file_name = outfile,
    pvOutputThreshold = p_val,
    useModel = model,
    errorCovariance = numeric(),
    verbose = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = TRUE,
    pvalue.hist = FALSE
  )
  message("Done: ", quest_qtl$param$output_file_name, "\n")
}

dir.create(out_mapping, recursive = TRUE, showWarnings = FALSE)

# start
chr <- sub("\\.txt$", "", basename(snp_path))   # e.g. chr1
message("### Now processing ", chr, " ... ###")

# Load covariates
message("Load covariate ", cov_path)
cov_all <- fread(cov_path, header = TRUE)
cov_all <- cov_all %>% column_to_rownames(names(cov_all)[1])

# split covariates and check existance in cov_all
cov_main <- str_split(cov_main_unsplit, ",")[[1]]
#cov_interact <- str_split(cov_interact_unsplit, ",")[[1]]
diff <- setdiff(cov_main, row.names(cov_all))
if (length(diff) != 0) {
  stop("Main covariate(s) ", paste(diff, collapse = ", "),
  " not found in covariate file")
}
if (cov_interact != "main") {
  diff <- setdiff(cov_interact, row.names(cov_all))
  if (length(diff) != 0) {
    stop("Interact covariate(s) ", paste(diff, collapse = ", "),
    " not found in covariate file")
  }
}

# Load genotype
message("Load genotype ", snp_path)
lb_all <- fread(snp_path, header = TRUE) %>% column_to_rownames(var = "SNP")
#names(lb_all) <- gsub("\\b\\d+_", "", names(lb_all))  # remove leading numbers before underscore with gsub
#names(lb_all) <- gsub("\\b\\d+_", "", names(lb_all))  # twice needed sometimes

# Load phenotype (already normalised)
message("Load phenotype ", pheno_path)
pheno_df <- fread(pheno_path, header = TRUE)
pheno_df <- pheno_df %>% column_to_rownames(names(pheno_df)[1])


info <- list()
message("Total ", pheno_type, ", : ", nrow(pheno_df))
i <- 0
pheno_name <- "pbmc_7d_il17_calbhy"
for (pheno_name in rownames(pheno_df)) {
  i <- i + 1
  pheno <- pheno_df %>% filter(rownames(.) == pheno_name)
  orig_len <- ncol(pheno)
  message("\nNow processing (", i, "/", nrow(pheno_df), "): ", pheno_name)
  outfile <- paste0(out_mapping, "/", pheno_name, ".tsv")

  pheno <- pheno %>%
    select_if(~ !any(is.na(.))) %>%  #remove without pheno measurement
    select(any_of(names(lb_all))) %>%  # select only samples that are genotyped
    select(any_of(names(cov_all)))  # only samples that have covariate info
  message("Sample change after matching sample names: ", orig_len, " -> ", ncol(pheno), "\n")

  # checking if this is chr 1 and main effect to only create this list once
  if (chr == "chr1" && cov_interact == "main") {
    info[pheno_name] <- ncol(pheno)
  }

  # check sample number, these samples will error at mapping, so skip
  if (ncol(pheno) <= (nrow(cov_all) + 1 + 1)) {
    message("Warning: ", pheno_name, " contains only ", ncol(pheno), " samples!!!! skipping...")
    file.create(outfile)
  } else {
    # select only samples that exist in this specific phenotype
    lb <- lb_all %>% subset(select = names(pheno))
    cov <- cov_all %>% subset(select = names(pheno))

    # mapping without interaction
    if (cov_interact == "main") { # aka no interaction
      cov_only_main <- cov %>% filter(rownames(.) %in% cov_main)
      cov_rows <- apply(cov_only_main, 1, function(row) length(unique(row)) == 1)
      if (any(cov_rows)) {  
        # if any covariates only contain one single value, cannot be mapped, still creating file for snakemake
        message("Warning: ", pheno_name, " has same covariate values for all samples, running without the following covariate(s)...")
        print(which(cov_rows))
        # exclude this covariate from cov_only_main
        cov_only_main <- cov_only_main %>% slice(-which(cov_rows))
        #file.create(outfile)
      } #else {
      quest_qtl <- perform_mapping(modelLINEAR, pheno_name, pheno, lb, cov_only_main, outfile)
      #}
    } else { 
      # make a new covariate df that includes only the 'main' covariates
      # and the interact covariate as last row
      cov_this <- cov %>%
        # filter on only main
        filter(rownames(.) %in% c(cov_main, cov_interact)) %>%
        # extra remove incase interact = main too
        slice(-which(row.names(.) == cov_interact)) %>%
        # then add the interact as last cov
        bind_rows(filter(cov, row.names(cov) == cov_interact)) %>%
        # exclude samples that do not have this covariate
        select_if(~ !any(is.na(.)))

      #pheno <- pheno %>% subset(select = names(cov_this))  # extra filter?? TODO
      cov_rows <- apply(cov_this, 1, function(row) length(unique(row)) == 1)
      if (any(cov_rows)) {  
        # if any covariates only contain one single value, cannot be mapped, still creating file for snakemake
        message("Warning: ", pheno_name, " has same covariate values for all samples, skipping...")
        file.create(outfile)
      } else {
        # also exclude these samples from lb and pheno
        lb_this <- lb %>% subset(select = names(cov_this))
        pheno_this <- pheno %>% subset(select = names(cov_this))

        # perform mapping with last covariate as interact
        quest_qtl <- perform_mapping(modelLINEAR_CROSS, pheno_name, pheno_this, lb_this, cov_this, outfile)
      }
    }  
  }
}

# idk why but this part seems to error randomly despite actually working so im adding this as failsafe..
tryCatch({
  if (length(info) > 0) {
    info_df <- stack(info)
    names(info_df) <- c("samples", pheno_type)
    fwrite(info_df, file = paste0(out_mapping, "/mapped_samples.txt"), sep = "\t")
    message("Saved sample info in ", paste0(out_mapping, "/mapped_samples.txt"))
  }

  message("Done with ", chr)
}, error = function(e) {
  message("idek ", e)
})
