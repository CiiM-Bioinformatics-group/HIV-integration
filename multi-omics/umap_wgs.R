suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
})

covs = '/vol/projects/CIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "/vol/projects/CIIM/2000HIV/cQTL/mofa/out/"

wgsf <- '/vol/projects/zzhang/projects/wp_2000hiv/inputs/phenotypes/230221_2000hiv_study_export_processed_EC.donor_with_wgs.csv'
wgs <- fread(wgsf)%>%dplyr::rename('sample' = ID)

model <- readRDS(paste0(outdir,'model_corrected_scaled.rds'))

covdf <- left_join(samples_metadata(model)%>%dplyr::select(sample), wgs, 'sample')
samples_metadata(model) <- covdf

df <- left_join(samples_metadata(model),model@dim_red$UMAP, 'sample')
ggplot(df,aes(x=UMAP1,y=UMAP2,color=EC_combined))+
  geom_point()+theme_bw()+theme(panel.grid = element_blank())
ggsave(paste0(outdir, 'umap_wgs',append,'.pdf'),width = 6,height = 4)

#Test factors that are different in this group
df_wilcox <- df %>% 
  mutate(wgs = ifelse(sample %in% wgs$sample, 'wgs','no'))%>%
  dplyr::select(c(wgs,colnames(.)[grepl('single_group.Factor',colnames(.))]))%>%
  pivot_longer(colnames(.)[grepl('single_group.Factor',colnames(.))], names_to = 'factor')%>%
  group_by(factor)%>%
  rstatix::wilcox_test(value ~ wgs)%>%
  mutate(padj = p.adjust(p, method = 'fdr'))
fwrite(df_wilcox, paste0(outdir, 'factor_wgs_wilcox',append,'.csv'))
