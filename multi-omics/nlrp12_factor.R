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
  library(rstatix)
})

covs = '/vol/projects/CIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "/vol/projects/CIIM/2000HIV/cQTL/mofa/out/"

append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Reading each correlation
f6_cor <- fread(cmd = paste0('grep -e Factor6 -e p.adj ',outdir,'factors_cov_cor',append,'.csv'))%>%melt()%>%
  mutate(variable = gsub('p.adj','padj', variable) %>% sub('\\.', ',', .))%>%
  separate(variable, into = c('var','variable'),sep = ',')%>%
  pivot_wider(names_from = 'var', values_from = 'value')
f6_clin <- fread(cmd = paste0('grep -e Factor6 -e p.adj ',outdir,'clin_cor/clin_df_cor',append,'.csv'))%>%melt()%>%
  mutate(variable = gsub('p.adj','padj', variable) %>% sub('\\.', ',', .))%>%
  separate(variable, into = c('var','variable'),sep = ',')%>%
  pivot_wider(names_from = 'var', values_from = 'value')%>%
  na.omit()%>%
  mutate('fdr' = p.adjust(p, method = 'fdr'))


#Basically correlate factor 6 to NLRP12
f6 <- get_factors(model, as.data.frame = T)
gex <- get_data(model, views = 'gex', as.data.frame = T)
prot <- get_data(model, views = 'prot', as.data.frame = T)

f6_nlrp12  <- inner_join(
  filter(f6, factor == 'Factor6'),
  filter(rbind(gex,prot), feature %in% c('NLRP12','NLRP3','CASP1','PARP1','IL1B','IL18','IL16','IL1RN','NFKB1',
                                         'IL1B_Inflammation', 'CASP1_Neurology','PARP1_Inflammation',
                                         'NFKB1_Cardiometabolic_II',
                                         'IL18_Inflammation','IL16_Inflammation','IL1RN_Inflammation')),
  by='sample'
)

f6_prots  <- inner_join(
  filter(f6, factor == 'Factor6'),
  prot,
  by='sample'
)%>%group_by(feature)%>%cor_test(value.x, value.y)%>%mutate(padj=p.adjust(p, method = 'fdr'))

ggplot(f6_nlrp12 %>% filter(grepl('^IL',feature))%>%
         mutate(feature = gsub('_Inflammation','',feature) %>% factor(levels = c('IL1B','IL1RN','IL18','IL6'))),
       aes(x=value.x,y=value.y))+
  geom_point(alpha=.1)+geom_smooth(method = 'lm')+
  facet_grid(view~feature)+
  stat_cor(aes(label =paste(..r.label..,
                            cut(..p.., breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")),
                            sep = "~")))+
           xlab('Factor6')+ylab('Scaled value')+
  theme_bw()+
  theme(panel.grid = element_blank(),strip.background = element_blank())
ggsave(paste0(outdir, 'dotplot_factor6_inflammasome_il_prot_gex.pdf'),width =5,height = 3)

ggplot(f6_nlrp12 %>% filter(!grepl('^IL',feature))%>%separate(feature,'feature'),
       aes(x=value.x,y=value.y))+
  geom_point(alpha=.1)+geom_smooth(method = 'lm')+
  facet_grid(view~feature)+
  stat_cor(aes(label =paste(..r.label..,
                            cut(..p.., breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")),
                            sep = "~")))+
  xlab('Factor6')+ylab('Scaled value')+
  theme_bw()+
  theme(panel.grid = element_blank(),strip.background = element_blank())
ggsave(paste0(outdir, 'dotplot_factor6_inflammasome_gex.pdf'),width =25/4,height = 3)


f6_nlrp12 %>% filter(grepl('IL1B',feature))%>%
  mutate(feature = gsub('_Inflammation','',feature))%>%
  pivot_wider(names_from = view, values_from = value.y)%>%
  mutate(IL1B_diff = prot-gex)%>%
  ggplot(aes(x=value.x, y=IL1B_diff))+
  geom_point()+geom_smooth(method = 'lm')+
  stat_cor()+
  xlab('Factor6')+
  theme_bw()+
  theme(panel.grid = element_blank())

#Get all QTLs for the top snp
allqtls <- fread(cmd =
                   'grep rs3974831 /vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/*/mapping/main_genomewide.tsv',
                 col.names = c('file','feature','beta','t','p'))%>%
  mutate(file = gsub('/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/','',file))%>%
  separate(file, into = 'view', sep = '/')
