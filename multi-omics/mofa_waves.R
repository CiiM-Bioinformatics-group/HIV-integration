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
  library(parallel)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(minfi)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"

append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

factors_df <- get_factors(model, as.data.frame = T)
ggplot(factors_df, aes(x=value+as.numeric(gsub('Factor','',factor))*7,fill=factor))+
  geom_density(alpha = .2,color='white')+
  #scale_fill_manual(values=c(rep(c('#E7B4C7','#DB8DAA'),10),'#E7B4C7'))+
  #facet_grid(~factor, switch = 'x', scales = 'free', space = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        panel.spacing = unit(0, "lines"), strip.background = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm'),strip.placement = 'outside', legend.position = 0)

