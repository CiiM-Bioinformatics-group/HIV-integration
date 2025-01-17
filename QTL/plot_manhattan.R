#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
library(ggman)
library(ggrastr)
library(CMplot)

# if package ‘ggman’ is not available for this version of R
#install.packages("remotes")
#remotes::install_github("mkanai/ggman")

args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
cov <- args[2]
genome_build <- args[3]
meta <- as.logical(args[4])
output <- args[5]
pval <- as.numeric(args[6])
if (meta) {
  input <- args[7]
} else {
  input <- as.vector(args)[-c(1:6)] # everything except the args that come before, dont put arguments after this
  input <- paste(input, collapse = " ")
}

# see if args are correct
message("path: ", path,
        "\ncov: ", cov,
        "\ngenome_build ", genome_build,
        "\nmeta: ", meta,
        "\noutput: ", output,
        "\npval: ", pval,
        "\ninput: ", input)

# Load in data, only taking p-value <0.05
gw <- fread(cmd = paste0("awk 'NR==1{print} NR>1{if($5 < 0.05) print}' ", input))

if (nrow(gw) > 0) {
  gw %>% arrange(`p-value`) %>% head

  # Add chr pos
  gw <- gw %>%
    rename( "P" = `p-value`) %>%
    group_by(SNP) %>%
    summarise(`P` = min(`P`)) %>%
    # remove row where first column is SNP
    filter(!grepl("SNP", SNP)) %>%
    separate(SNP, into = c("CHR", "BP"), sep = ":",
      remove = FALSE, convert = TRUE, extra = "drop") %>%
    mutate(CHR = gsub("chr", "", CHR)) # remove the chr in chromosome if it's there

  #remove chromosome 23 as we are not interested in it currently
  gw <- gw[gw$CHR != "X", ]
  gw <- gw[gw$CHR != "23", ]
  #gw$CHR <- ifelse(gw$CHR == "X", 23, gw$CHR) # change X to 23
  #gw$CHR <- as.numeric(gw$CHR) # only after changing X to 23 otherwise it'll be NA

  outdir <- path
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)  # Create plot dir in case it does not exist yet
  setwd(outdir)

  tryCatch({
    CMplot::CMplot(gw, type = "p", plot.type = "m", band = 0.5, LOG10 = TRUE,
      threshold = c(5e-8, pval),
      threshold.lty = 2, threshold.lwd = 1, threshold.col = "red",
      amplify = TRUE, width = 10, height = 5,
      signal.col = NULL, chr.den.col = NULL, file = "pdf", file.output = TRUE,  # memo = cov,
      verbose = TRUE, cex = 0.8, main = cov, dpi = 200,
      #col = colors
      ) # “jpg”, “pdf”, “tiff”

    # use system command to change output file name
    # might cause conflicts if multiple manhattan plots are made at same time
    message(paste0("mv ", outdir, "/Rect_Manhtn.P.pdf ", output))
    system(paste0("mv ", outdir, "/Rect_Manhtn.P.pdf ", output))

    man_gg <- topr::manhattan(gw %>% dplyr::select(CHR,BP,P)%>%setNames(c('CHROM','POS','P')),
                      annotate = pval, sign_thresh = pval)+ylab('')+
        xlab('')+
        theme(axis.text.x = element_blank(), axis.line.x = element_blank(),axis.ticks.x = element_blank(),
              plot.margin = unit(c(0,0,0,0),'mm'))

    ggsave(gsub('.pdf','.png',output), man_gg, width = 10, height = 5, dpi = 'print')


  }, error = function(e) {
    message("CMplot failed: ", e)
  })
} else {
  message("No significant results")
  # make empty plot
  pdf(output, width = 10, height = 5)
}
