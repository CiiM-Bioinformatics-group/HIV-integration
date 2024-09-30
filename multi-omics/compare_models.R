suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
  library(reticulate)
  library(ggExtra)
  library(patchwork)
  library(ComplexHeatmap)
  library(clusterProfiler)
})

outdir = "2000HIV/cQTL/mofa/out/"

model_corrected_scaled <- readRDS(paste0(outdir,"model_corrected_scaled.rds"))
model_uncorrected_scaled <- readRDS(paste0(outdir,"model_uncorrected_scaled.rds"))
model_corrected_unscaled <- readRDS(paste0(outdir,"model_corrected_unscaled.rds"))
model_uncorrected_unscaled <- readRDS(paste0(outdir,"model_uncorrected_unscaled.rds"))

all_models <- list('both' = model_corrected_scaled, 
                   'scaled' = model_uncorrected_scaled,
                  'corrected'  = model_corrected_unscaled,
                  'neither' = model_uncorrected_unscaled)

compare_factors(all_models, show_colnames=F, treeheight_row=0, treeheight_col=0)
compare_factors(all_models[c(1,3)], show_colnames=F, treeheight_row=0, treeheight_col=0)


#Variances explained
variance_barplot <- function(mofa){
  variance_perFactor = mofa@cache$variance_explained$r2_per_factor
  variance_perFactorMelted = melt(variance_perFactor)
  ggplot(variance_perFactorMelted, aes(x=Var1, y=value,fill=Var2))+geom_col(position = 'stack')+theme_bw()+
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90))+ylab('Variance explained (%)')+xlab('Factors')+
    scale_fill_manual(values = viewcol)
}

variancesplot <- lapply(all_models,variance_barplot)
