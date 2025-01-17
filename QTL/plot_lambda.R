#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
outfile <- args[1]
mapping_out <- args[2]
input <- as.list(args)[-c(1:3)] # everything except the other args

message("outfile: ", outfile, 
        "\nmapping_out: ", mapping_out,
        "\ninput: ", paste(input, collapse = ", "))

# create outdir
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

# read values from files
lambda <- do.call(rbind, lapply(input, function(file) {
    gene <- tools::file_path_sans_ext(basename(file))
    lam <- as.numeric(readLines(file))
    data.frame(gene = gene, lam = lam)
}))

# plot
set.seed(123)
lambda_plot <- ggplot(lambda, aes(x = "gene", y = lam)) +
    geom_boxplot() + geom_jitter() + theme_bw()
# add name to plot if lambda >1.1 or <0.9
if (nrow(lambda[lambda$lam > 1.1 | lambda$lam < 0.9, ]) > 0) {
    lambda_plot <- lambda_plot + geom_text(data = lambda[lambda$lam > 1.1 | lambda$lam < 0.9, ],
    aes(label = gene), position = "jitter")
}

# save plot
ggsave(plot = lambda_plot, filename = outfile, width = 6, height = 10)
message(outfile)

# function for making qqplot
gq <- function(lambda, type, log_p = TRUE) {
    # read in summary statistics for these "genes"
    summary <- fread(cmd = paste0("awk 'NR == 1 || FNR > 1' ", mapping_out, "/chr", 1:22, "/", lambda$gene, ".tsv"))
    names(summary) <- c("SNP", "gene", "beta", "t-stat", "p")

    x <- summary$p
    x <- x[!is.na(x) & !is.nan(x) & !is.null(x) & is.finite(x)]
    if (log_p) {
        x <- -log10(x[x < 1 & x > 0])
    }

    data <- data.frame(observed = sort(x, decreasing = TRUE),
                        expected = -log10(ppoints(length(x))))
    plt <- ggplot(data, aes(expected, observed)) +
        geom_point(size = 0.02) +
        geom_abline(slope = 1, colour = "red") +
        xlab(expression(Expected ~ ~ -log[10](italic(P)))) +
        ylab(expression(Observed ~ ~ -log[10](italic(P)))) +
        ggtitle(paste0("qqplot ", type, " lambda = ", lambda$lam))

    outname <- paste0(dirname(outfile), "/", type, "_", lambda$gene, ".png")
    ggsave(plot = plt, filename = outname, width = 4, height = 4)
    message(outname)
}

# get min & max lambda, and plot qqplots for these
min_lambda <- lambda %>% filter(lam == min(lambda$lam))
gq(min_lambda, "min")

max_lambda <- lambda %>% filter(lam == max(lambda$lam))
gq(max_lambda, "max")
