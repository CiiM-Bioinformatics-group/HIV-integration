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
  library(patchwork)
  library(corrplot)
  library(rstatix)
})

#Inflammation score
covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'figure1'))
append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Colors: From jbgb13/peRReo
source('2000HIV/cQTL/mofa/code/peRReo.R')
viewcol <- latin_palette('buenavista', n=5)
names(viewcol) <- views_names(model)


prot <- get_data(model, views = 'prot', as.data.frame = T)
gex <- get_data(model, views = 'gex', as.data.frame = T)
cyt <- get_data(model, views = 'cyt', as.data.frame = T)

#Calculate a per-individual inflammatory score
#Overall inflammatory response to pathogens

#BTMs
btm <- read.gmt(paste0(outdir,'BTM_for_GSEA_20131008.gmt'))
inflam_genes <- btm %>% filter(grepl('\\(M29\\)', term)) %>% pull(gene)

#Lets go for three cytokines
inflam_genes <- c('IL1B')

#Based on gen names
prot_inflam <- prot[grepl(paste(paste0('^',inflam_genes,'_'), collapse = '|'), prot$feature),]%>%
  separate(feature, into = c('feature'), sep = '_') %>% filter(feature %in% inflam_genes) %>%
  group_by(sample) %>% summarise('prot_score'=  mean(value))

#Based on averaging all of 'inflammation' genes
#prot_inflam <- prot[grepl('Inflammation', prot$feature),]%>%
#  group_by(sample) %>% summarise('prot_score'=  mean(value))


#Same stuff for gex
gex_inflam <- gex %>% filter(feature %in% inflam_genes)%>%
  group_by(sample) %>% summarise('gex_score'=  mean(value))

#For cyt
cyt_inflam <- cyt[grepl('il1b', cyt$feature),]%>%
  group_by(sample) %>% summarise('cyt_score'=  mean(value))

cyt_all <- cyt[grepl('il1b', cyt$feature),]

#gex divided by prot
gex_prot_inflam <- inner_join(gex_inflam,prot_inflam)%>%mutate('gex_prot_score' = gex_score-prot_score)

all_inflam <- inner_join(gex_prot_inflam,cyt_inflam)%>%
  column_to_rownames('sample')%>%
  mutate_all(scale)

all_inflam %>%
  pivot_longer(colnames(.), names_to = 'score')%>%
  mutate(score=gsub('_score','',score)%>%factor(levels = c('cyt','prot','gex')))%>%
  filter(score %in% names(viewcol))%>%
  ggplot(aes(x=value,fill=score))+
  geom_density(color='white')+
  facet_wrap(~score, ncol = 1)+
  scale_fill_manual(values = viewcol)+
  xlab('score')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90, hjust=1,vjust=0.2),
        panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), panel.border = element_blank(),axis.line.x = element_line(),
        plot.margin = unit(c(0,0,0,0), 'cm'),strip.placement = 'outside', legend.position = 0)

ggsave(paste0(outdir, 'figure1/il1b_score_density.pdf'),width=3,height=4)

ggplot(all_inflam, aes(x=prot_score,y=cyt_score))+
  geom_point()+geom_smooth(method = 'lm')+
  theme_bw()+theme(panel.grid = element_blank())+
  stat_cor()
ggsave(paste0(outdir, 'figure1/il1b_score_prot_cyt_cor.pdf'),width=3,height=3)

ggplot(all_inflam, aes(x=gex_score,y=cyt_score))+
  geom_point()+geom_smooth(method = 'lm')+
  theme_bw()+theme(panel.grid = element_blank())+
  stat_cor()
ggsave(paste0(outdir, 'figure1/il1b_score_gex_cyt_cor.pdf'),width=3,height=3)


ggplot(all_inflam, aes(x=gex_score,y=prot_score))+
  geom_point()+geom_smooth(method = 'lm')+
  theme_bw()+theme(panel.grid = element_blank())+
  stat_cor()
ggsave(paste0(outdir, 'figure1/il1b_score_gex_prot_cor.pdf'),width=3,height=3)

ggplot(all_inflam, aes(x=gex_prot_score,y=cyt_score))+
  geom_point()+geom_smooth(method = 'lm')+
  theme_bw()+theme(panel.grid = element_blank())+
  stat_cor()
ggsave(paste0(outdir, 'figure1/il1b_score_gex_prot_cyt_cor.pdf'),width=3,height=3)

#CV for the different layers
CV <- function(x){
  if(any(x<0)) x <- x + abs(min(x,na.rm = T))
  sd(x, na.rm = T)/mean(x, na.rm=T)
}

prot_cv <- CV(prot_inflam$prot_score)
gex_cv <- CV(gex_inflam$gex_score)
cyt_cv <- cyt_all %>% group_by(feature) %>% summarise('CV' = CV(value)) %>% pull(CV, name = feature)


#One plot with all correlations and variation
variance_perFactorMelted = reshape2::melt(model@cache$variance_explained$r2_per_factor)

barplot_filled <- variance_perFactorMelted %>%
  mutate(Var2 = factor(Var2, levels = c('meth','gex','prot','metab','cyt')),
         Var1 = factor(Var1, levels = paste0('Factor',21:1)))%>%
  ggplot(aes(y=Var1, x=value,fill=Var2))+geom_col(position = 'fill')+theme_bw()+
  scale_y_discrete(position = 'right')+
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),axis.text.y.right = element_text(hjust = 0),
        axis.title = element_blank(), axis.text.x = element_blank(),panel.border = element_blank(),strip.placement = 'outside',
        panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm'),legend.title = element_blank())+
  scale_fill_manual(values = viewcol)

barplot_black <- variance_perFactorMelted %>%
  group_by(Var1)%>%summarise(value=sum(value))%>%
  mutate(Var1 = factor(Var1, levels = paste0('Factor',21:1)))%>%
  ggplot(aes(y=Var1, x=value))+geom_col(fill='darkgrey')+theme_bw()+
  xlab('')+
  facet_wrap(~'')+
  scale_x_reverse(position = 'top')+
  scale_y_discrete(position = 'right')+
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),axis.text.y.right = element_text(hjust = 0),
        axis.title = element_blank(), panel.border = element_blank(),strip.placement = 'outside',
        panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm'),legend.title = element_blank(),
        axis.line.x.top = element_line())+
  scale_fill_manual(values = viewcol)


#Reading the cors

corcsv_to_plot <- function(fulltable){
  #fulltable <- fread(file)%>%column_to_rownames('V1')

  #Those with one significant hit
  onesig <- fulltable %>% group_by(clin_var) %>% summarise('sig' = any(p<0.05))%>%
    filter(sig)%>%pull(clin_var)

  fulltable <- fulltable %>% filter(clin_var %in% onesig)%>%
    mutate(sign = ifelse(is.na(estimate), 1, sign(estimate)))%>%
    mutate(value = -log10(padj) * sign)%>%
    mutate(factor = factor(factor, levels = paste0('Factor',1:length(unique(fulltable$factor)))))


  ggplot()+
    geom_point(data=fulltable%>%
                 mutate(size=ifelse(abs(value) > -log10(0.001), '4', ifelse(abs(value)>-log10(0.01),'3',
                                                                               ifelse(abs(value)>-log10(0.05),'2',
                                                                                      ifelse(p<0.05,'1','0'))))),
               aes(y=factor,x=clin_var,color=value,size=size))+
    geom_text(data = fulltable%>%
                mutate(value=ifelse(abs(value) > -log10(0.001), '***', ifelse(abs(value)>-log10(0.01),'**',
                                                                 ifelse(abs(value)>-log10(0.05),'*',
                                                                        ifelse(p<0.05,' ',''))))),
              aes(y=factor,x=clin_var,label=value),vjust=.75)+
    scale_color_gradient2( low = 'darkblue',high = 'darkred', name = '-log10(P)')+
    scale_size_discrete(name = 'Significance')+
    facet_grid(~type,scales = 'free', space = 'free')+
    scale_y_discrete(limits=rev)+
    scale_x_discrete(position = 'top')+
    theme_bw()+theme(panel.grid = element_blank(),panel.border = element_blank(),
                     axis.ticks = element_blank(),axis.text.x = element_text(angle=45, hjust = 0),
                     axis.title = element_blank(),strip.placement = 'outside',
                     panel.spacing = unit(.2, "cm"), strip.background = element_blank(),
                     plot.margin = unit(c(0,0,0,0), 'cm'))+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=.5)

}

factors_all_test <-  rbind(
  fread(paste0(outdir,'cov_cor/factors_cov_test',append,'.csv'))%>%mutate('type' = 'Covariates'),
  fread(paste0(outdir,'clin_cor/clin_df_test',append,'.csv'))%>%
    separate(clin_var, into = c('type','clin_var'),sep = '%')%>%
    mutate(clin_var = ifelse(is.na(clin_var), type, clin_var))%>%
    mutate(type = ifelse(type == 'PLAQUE_YESNO', 'CV', type))
)%>% mutate(type = ifelse(clin_var %in% c('INR_YESNO','RP_YESNO','CD4_LATEST', 'CD8_LATEST', 'CD4_NADIR', 'HIV_DURATION'),
                          'HIV', type))

factors_all_test_main <- factors_all_test %>%
  subset(!clin_var %in%  c('`CD4T/CD8T`', '`CD4_LATEST/CD8_LATEST`', 'CD4T', 'CD8T',
                          'cplvvaccinated', 'ETHNICITY.Asian','SMOKING_PY','Not_done'))%>%
  subset(!grepl('AGE|BMI|CANNABIS|CMV|Institute|ETHNICITY|EC_|SEX_BIRTH',clin_var))%>%
  mutate(clin_var = gsub('_YESNO','',clin_var)%>%gsub('MH_TR_','',.)%>%gsub('_',' ', .)%>%str_to_sentence(),
         type = gsub('MH_TR_','',type)%>%gsub('_D','', .)%>%str_to_sentence())


factor_all_plot <- corcsv_to_plot(factors_all_test_main)
clin_vars_main <- factor_all_plot[["layers"]][[1]][["data"]][["clin_var"]] %>% unique()


all_inflam <- dplyr::select(all_inflam,-c(gex_prot_score))

#Factors that capture this general IL-1B inflammation? Probably age but what else
factor_df <- get_factors(model, as.data.frame = T)
factor_il1b <- inner_join(factor_df, all_inflam %>%
                            rownames_to_column('sample')%>%
                            pivot_longer(colnames(.)[-1], names_to = 'score')%>%
                            mutate(score=gsub('_score','',score)%>%factor(levels = c('cyt','prot','gex'))),
                          by = 'sample')

factor_il1b_cor <- factor_il1b %>%
  group_by(factor,score)%>%
  cor_test(value.x,value.y)%>%
  mutate(padj = p.adjust(p, method = 'fdr'))

fwrite(factor_il1b_cor, paste0(outdir,'figure1/factors_to_il1b_cor',append,'.csv'))

factor_il1b_cor_hm <- factor_il1b_cor %>%
  mutate(cor = ifelse(padj<0.05,cor, 0))%>%
  mutate(factor=factor(factor, levels = paste0('Factor',21:1)))%>%
  ggplot(aes(y=factor,x=score,fill=cor))+
  geom_tile()+
  scale_fill_gradient2(limits=c(-1,1), low = 'darkblue',high = 'darkred')+
  xlab('')+
  scale_x_discrete(position = 'top')+
  facet_wrap(~'')+
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),axis.text.y = element_blank(),
        axis.text.x.top = element_text(angle = 90, vjust=0),
        axis.title.y = element_blank(), panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm'),legend.position = 0,strip.placement = 'outside',
        panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        panel.background = element_blank())


barplot_black+theme(legend.position = 0,axis.text.y.right = element_blank())+
  barplot_filled+theme(legend.position = 0,axis.text.y = element_blank())+
  #factor_eaa_cor_hm+
  factor_il1b_cor_hm+
  factor_all_plot+theme(axis.text.y = element_blank(),legend.position = 0)+
  plot_layout(nrow=1, widths = c(5,5,
    #1,
    3,5*2+5*9.6))
ggsave(paste0(outdir,'figure1/variance_cov_clin_cor_superplot',append,'.pdf'), width = 12.5, height = 6)

#Legends
barplot_filled_legend <- as_ggplot(get_legend(
  barplot_filled+theme(legend.position = 'top')+guides(fill = guide_legend(reverse = TRUE))
  ))
ggsave(plot = barplot_filled_legend,paste0(outdir,'figure1/variance_cov_clin_cor_superplot_barplot_legend',append,'.pdf'), width = 3.5, height = 1)



#The genetic component
mofaqtl <- fread('meta_cQTL/out/2000HIV-EU-discovery/mofa/mapping/main_1e-5.tsv')%>%
  group_by(gene) %>%slice_min(`p-value`,n = 1)

barplot_mofaqtl <- mofaqtl %>%
  mutate(gene = factor(gene, levels = paste0('Factor',21:1)))%>%
  ggplot(aes(y=gene, x=-log10(`p-value`),fill=-log10(`p-value`)))+
  geom_col()+theme_bw()+
  scale_fill_gradient(low = 'ivory2', high = '#6A4D74')+
  geom_vline(aes(xintercept=-log10(5e-8/21)), linetype = 'dashed')+
  xlab('')+
  facet_wrap(~'')+
  scale_x_continuous(position = 'top')+
  scale_y_discrete(position = 'right')+
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),axis.text.y = element_blank(),
        axis.title.y = element_blank(), panel.border = element_blank(),legend.position = 0,strip.placement = 'outside',
        panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        axis.line.x = element_line(),plot.margin = unit(c(0,0,0,0), 'cm'))

barplot_black+theme(legend.position = 0,axis.text.y.right = element_blank())+
  barplot_filled+theme(legend.position = 0)+
  #factor_eaa_cor_hm+
  factor_il1b_cor_hm+
  barplot_mofaqtl+factor_all_plot+theme(legend.position = 0,axis.text.y.left = element_blank())+
  plot_layout(nrow=1, widths = c(5,5,
    #1,
    3,5*2,5*9.6))
ggsave(paste0(outdir,'figure1/variance_cov_clin_cor_genetics_superplot',append,'.pdf'), width = 12.5, height = 6)


#Have we found a pathway that is different depending on the ART treatment?
#ifn_genes <- btm %>% filter(grepl('type I interferon response', term)) %>% pull(gene)
ifn_genes <- btm %>% filter(grepl('antiviral IFN', term)) %>% pull(gene)
#ifn_genes <- btm %>% filter(grepl('viral sensing \\& immunity\\; IRF2', term)) %>% pull(gene) %>% unique()

#Based on gen names
prot_ifn <- prot[grepl(paste(paste0('^',ifn_genes,'_'), collapse = '|'), prot$feature),]%>%
  separate(feature, into = c('feature'), sep = '_') %>% filter(feature %in% ifn_genes) %>%
  group_by(sample) %>% summarise('prot_score'=  mean(value))

#Same stuff for gex
gex_ifn <- gex %>% filter(feature %in% ifn_genes)%>%
  group_by(sample) %>% summarise('gex_score'=  mean(value))

all_ifn <- left_join(gex_ifn,prot_ifn)%>%
  column_to_rownames('sample')%>%
  mutate_all(scale)
ggplot(all_ifn, aes(x=prot_score,y=gex_score))+geom_point()+stat_cor()+geom_smooth()+
  theme_bw()+theme(panel.grid = element_blank())


#Factor 8 to all the disease endpoints
factor8_test <- factors_all_test %>% filter(factor == 'Factor8') %>%
  filter(grepl('MH_TR',type) | grepl('MH_TR',clin_var)) %>%
  mutate(type = ifelse(type == '', clin_var,type))%>%
  mutate(type = gsub('MH_TR_','',type)%>% gsub('_D','',.))%>%
  mutate(clin_var = gsub('MH_TR_','',clin_var))%>%
  mutate(clin_var = ifelse(clin_var == 'Other',paste0('Other_',type),clin_var))%>%
  mutate(range = conf.high-conf.low)

ggplot(factor8_test,aes(y=reorder(clin_var,-p)))+
  geom_point(aes(x=estimate))+
  geom_linerange(aes(xmin=conf.low,xmax=conf.high))+
  geom_vline(aes(xintercept=0), linetype='dashed')+
  facet_grid(type~., space = 'free',scales = 'free_y')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank(), axis.lines = element_line())+
  ylab('Clinical variable')+xlab('Estimate')
ggsave(paste0(outdir,'figure1/Factor8_clinical',append,'.pdf'), width = 3, height = 5)
