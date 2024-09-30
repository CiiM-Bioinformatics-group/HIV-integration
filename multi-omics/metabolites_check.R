suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
})

#DECISION: Leave PC1 in the metabolites as it seems the cause is not clear and MOFA does find it and therefore does not
#affect our analyses

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"

append = '_corrected_scaled'

model <- readRDS(paste0(outdir,'model', append, '.rds'))

metab_unscaled <- fread('2000HIV/Metabolites/2000HIV_GeneralMetabolomics_raw_metabolome_data_outliers_and_pSS_excluded_1902samples_touse.csv')%>%
  column_to_rownames(colnames(.)[1]) %>% .[-c(1,2,3)]

metabs <- get_data(model, views = 'metab')[[1]][[1]]
metabs <- metabs[!duplicated(rownames(metabs)),]
covs <- samples_metadata(model)

metab_pc <- prcomp(t(metabs))$x %>% as.data.frame() %>% rownames_to_column('sample')
metab_pca <- prcomp(t(metabs))

factors <- get_factors(model)[[1]] %>% as.data.frame() %>% rownames_to_column('sample')
metabs_df <- metab_unscaled %>% as.data.frame()%>% rownames_to_column('sample')
#dups <-  colnames(metabs_df) %in% colnames(metabs_df)[duplicated(colnames(metabs_df))]
#colnames(metabs_df)[dups] <- paste0(colnames(metabs_df)[dups],'_',1:sum(dups))

df <- left_join(covs, metab_pc)%>%left_join(factors)%>%left_join(metabs_df)

ggplot(df, aes(x=PC1,y=PC2,color=days_post_lockdown))+geom_point()
ggplot(df, aes(x=PC1,y=PC2,color=cplv))+geom_point()
ggplot(df, aes(x=PC1,y=PC2,color=Factor1))+geom_point()+
  scale_color_gradient2()
cor.test(df$PC1, df$Factor1)

metab_loads <- prcomp(t(metabs))
metabs_pc_top <- abs(metab_loads[["rotation"]][,'PC1']) %>% sort(decreasing = T) %>%
  .[1:20] %>%
  names()
metab_ids  <- fread('2000HIV/Metabolites/ionMz2ID.tsv') %>%
  mutate(ionMz = gsub('X','',ionMz))%>% pull(ionMz, name=name)
metabs_pc_top_name <- unname(metab_ids[metabs_pc_top])

df_top_long <- df %>% dplyr::select(PC1, metabs_pc_top_name)%>%
  pivot_longer(metabs_pc_top_name)

ggplot(df_top_long, aes(x=PC1,y=log1p(value)))+geom_point()+
  facet_wrap(~name)+
  scale_color_gradient2(midpoint = 1, mid = 'white', low = 'blue')+
  geom_smooth()+stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())

df_long <- df %>% 
  pivot_longer(colnames(metabs_df)[-1])%>%
  group_by(name)%>%
  rstatix::cor_test(PC1, value)

#Associate factor1 to all covariates
df_pc1 <- df %>% dplyr::select(PC1, colnames(samples_metadata(model))[-c(1:2)]) %>%
  dplyr::select(-c('ETHNICITY.White','Vaccine_Group','cov','vac','iso','iso2','days_post_cov','days_post_vacc',
                   'ETHNICITY'))%>%
  mutate_at('CMV_IgG_IU.mL', as.numeric)#Remove things with 1 instance
lm_pc1 <- lm(PC1 ~ ., data = df_pc1)
summary(lm_pc1)

ggplot(df_pc1, aes(x=PC1,y=days_post_lockdown))+geom_point()+
  geom_smooth()+stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())

ggplot(df_pc1, aes(x=PC1,y=season_cos))+geom_point()+
  geom_smooth()+stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())

lm_pc1_2 <- lm(PC1 ~ season_cos+days_post_lockdown, data = df_pc1)
summary(lm_pc1_2)

ggplot(df, aes(x=PC1,y=Factor1))+geom_point()+
  geom_smooth()+stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())

#Shall we then correct PC1 in the metabolites? Potentially killing the variation in some
metabs_correct_pc <- prcomp(metab_unscaled, center = T, scale = T)$x %>% as.data.frame() %>%
  rownames_to_column('sample') %>% .[,c(1,2)]

