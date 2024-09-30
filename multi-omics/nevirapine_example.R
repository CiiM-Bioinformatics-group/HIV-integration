suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
  library(patchwork)
  library(corrplot)
  library(rstatix)
})

#Inflammation score

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"

append = '_corrected_scaled'

#Load
clin_df <- fread(paste0(outdir,'clin_cor/clin_df',append,'.csv'))%>%dplyr::rename('Record.Id' = colnames(.)[1])%>%
  left_join(fread(covs), by = 'Record.Id')%>%
  column_to_rownames(colnames(.)[1])
colnames(clin_df) <- gsub('%', '_', colnames(clin_df))

clin_df_validation <- fread(paste0(outdir,'validation/clin_df_validation_alpha0.01',append,'.csv'))%>%dplyr::rename('Record.Id' = colnames(.)[1])%>%
  left_join(fread(covs), by = 'Record.Id')%>%
  column_to_rownames(colnames(.)[1])
colnames(clin_df_validation) <- gsub('%', '_', colnames(clin_df_validation))


ggboxplot(clin_df,
          x='CART_NNRTI_Nevirapine', y='Factor15', add = 'jitter', add.params = list(alpha=.05))+
  theme(strip.background = element_blank())+
  stat_compare_means(label.y = 5)
ggsave(paste0(outdir,'figure1/discovery_boxplot_nevirapine_Factor15.pdf'),width = 3,height = 3)

ggboxplot(clin_df_validation,
          x='CART_NNRTI_Nevirapine', y='Factor15', add = 'jitter', add.params = list(alpha=.05))+
  theme(strip.background = element_blank())+
  stat_compare_means(label.y = 1.5)
ggsave(paste0(outdir,'figure1/validation_boxplot_nevirapine_Factor15.pdf'),width = 3,height = 3)

ggboxplot(clin_df%>%filter(!is.na(ETHNICITY.Black)),
          x='CART_NNRTI_Nevirapine', y='Factor15', add = 'jitter', add.params = list(alpha=.05))+
  facet_wrap(~ETHNICITY.Black,ncol = 1,strip.position = 'right')+
  theme(strip.background = element_blank())+
  stat_compare_means(label.y = 5)
ggsave(paste0(outdir,'figure1/discovery_boxplot_nevirapine_ethnicity_black_Factor15.pdf'),width = 3,height = 3)

ggboxplot(clin_df_validation%>%filter(!is.na(ETHNICITY.Black)),
          x='CART_NNRTI_Nevirapine', y='Factor15', add = 'jitter', add.params = list(alpha=.05))+
  facet_wrap(~ETHNICITY.Black,ncol = 1, strip.position = 'right')+
  stat_compare_means(label.y = 1.5)+theme(strip.background = element_blank())
ggsave(paste0(outdir,'figure1/validation_boxplot_nevirapine_ethnicity_black_Factor15.pdf'),width = 3,height = 3)

#Plot the nevirapine metabolite
nv_metab <- 'Hydroxynevirapine glucuronide'
model <- readRDS(paste0(outdir,'model', append, '.rds'))
metabs <- get_data(model, view = 'metab', as.data.frame=T) %>%
    filter(feature == nv_metab)

metabs_df <- left_join(metabs, clin_df %>% rownames_to_column('sample'))
head(metabs_df)
ggplot(metabs_df, aes(y=value, x=Factor15))+
    geom_point()+
    geom_smooth(method='lm')+
    stat_cor()+
    theme_bw()+
    theme(panel.grid=element_blank())+
    ylab(nv_metab)
ggsave(paste0(outdir,'figure1/discovery_nevirapine_metab_cor_Factor15.pdf'),width = 3,height = 3)

ggboxplot(metabs_df,y='value', x='CART_NNRTI_Nevirapine',
    add = 'jitter', add.params = list(alpha=.05))+
    stat_compare_means(label.y = 1.5)+
    ylab(nv_metab)
ggsave(paste0(outdir,'figure1/discovery_nevirapine_metab_boxplot.pdf'),width = 3,height = 3)


#Check Regim or Regimen
clin_df_art <- clin_df %>%
  .[,grepl('^CART|Factor15',colnames(.))] %>%
  .[,!grepl('CART_CCR5_inhibitor',colnames(.))]

#clin_df_art%>%.[-c(1:7,ncol(.))] %>% Heatmap()

nevi_regimen <- clin_df_art %>% dplyr::select(colnames(.)[grepl('Nevira|Emtri|TAF|Lami|Aba', colnames(.))], Factor15)%>%
  mutate('Regimen' = ifelse(CART_NRTI_Abacavir == 1, 'Aba', ''),
         'Regimen' = ifelse(CART_NRTI_Lamivudine == 1, paste0('Lami + ', Regimen), Regimen),
         'Regimen' = ifelse(CART_NTRTI_TAF == 1, paste0('TAF + ', Regimen), Regimen),
         'Regimen' = ifelse(CART_NRTI_Emtricitabine== 1, paste0('Emtri + ', Regimen), Regimen),
         'Regimen' = ifelse(CART_NNRTI_Nevirapine == 1, paste0('Nevi + ', Regimen), 'No-Nevi'),
         'Regimen' = gsub(' \\+ $', '', Regimen))%>%
  filter(Regimen != 'Nevi')%>%
  mutate(Regimen = gsub(' \\+ ', '\\+\n', Regimen))%>%
  mutate(Regimen = factor(Regimen, c('No-Nevi','Nevi+\nEmtri', 'Nevi+\nEmtri+\nTAF', 'Nevi+\nLami', 'Nevi+\nLami+\nAba')))

ggboxplot(nevi_regimen, x = 'Regimen', y = 'Factor15', add = 'jitter', add.params = list(alpha=0.05))+
  stat_pwc(label = 'p.signif', hide.ns = T)
ggsave(paste0(outdir,'figure1/discovery_nevirapine_Factor15_regimens_boxplot.pdf'),width = 5,height = 3)
