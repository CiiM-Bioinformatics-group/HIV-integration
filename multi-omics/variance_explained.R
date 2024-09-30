suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
  library(reticulate)
  library(ggExtra)
  library(patchwork)
  library(ComplexHeatmap)
  library(patchwork)
  library(corrplot)
  library(rstatix)
  library(caret)
})


covs = '/vol/projects/CIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "/vol/projects/CIIM/2000HIV/cQTL/mofa/out/"

append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Only eur
eur_samples <- samples_names(model)[[1]][model@samples_metadata$ETHNICITY.White %in% 1]
model <- subset_samples(model, eur_samples)

#Colors: From jbgb13/peRReo
source('/vol/projects/CIIM/2000HIV/cQTL/mofa/code/peRReo.R')
viewcol <- latin_palette('buenavista', n=5)
names(viewcol) <- views_names(model)

#Feature selection: 
#1. 1e-5 for SNPs, FDR<0.05 for covariates and all data layers
cyt <- get_data(model, views = 'cyt', as.data.frame = T)

cyt1 <- filter(cyt, feature == 'pbmc_24h_il1b_lps')

#TODO: Should correct for covariates before selecting, otherwise we get age-related factors
predictors <- lapply(c('prot','meth','metab','gex'), function(mol){
  sig_mol <- get_data(model, views = mol)[[1]][[1]][,cyt1$sample] %>%
    apply(1, function(x) cor.test(x,cyt1$value)$p.value)%>%
    p.adjust(method = 'fdr') %>% .[.<0.05] %>% sort() %>% names()
  
  cor_mol <- get_data(model, views = mol)[[1]][[1]][sig_mol,]%>%
    t() %>% cor()
  
  nonsig_mol <- findCorrelation(cor_mol, cutoff = .4, names = T,exact = T)
  
  sig_mol <- sig_mol[!sig_mol %in%  nonsig_mol]
  
  data <- get_data(model, views = mol)[[1]][[1]][sig_mol,cyt1$sample]
})
  
var <- calculate_variance_explained(model)

