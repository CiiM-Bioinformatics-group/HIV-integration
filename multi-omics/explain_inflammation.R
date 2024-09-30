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
  library(lavaan)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'figure2'))
append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Colors: From jbgb13/peRReo
source('2000HIV/cQTL/mofa/code/peRReo.R')
viewcol <- latin_palette('buenavista', n=5)
names(viewcol) <- views_names(model)


#Explained variance

#Discovery:
#Y
get_il1b <- function(df) df %>% .[grepl('IL1B|ENSG00000125538', rownames(.),ignore.case = T ),] %>% t()

model_unscaled <- load_model(paste0(outdir,'model_corrected_unscaled.hdf5'))

prot_Y <- get_data(model_unscaled, views = 'prot', add_intercept = F, denoise = F)[[1]][[1]] %>% get_il1b()
gex_Y <- get_data(model_unscaled, views = 'gex', add_intercept = F, denoise = F)[[1]][[1]] %>% get_il1b()
cyt_Y <- get_data(model_unscaled, views = 'cyt', add_intercept = F, denoise = F)[[1]][[1]] %>% t()

Y <- cbind(data.frame('IL1B_protein' = t(prot_Y), 'IL1B_gene' = t(gex_Y)),cyt_Y)

#Predictors
covar_X <- samples_metadata(model) %>% dplyr::select(sample,AGE,SEX_BIRTH,BMI_BASELINE,season_sin,season_cos) %>%
  column_to_rownames('sample')
factor_X <- get_factors(model)[[1]]

X <- cbind(covar_X,factor_X)

#Validation
#Y
read_validation <- function(view){
  fread(paste0(outdir,'validation/',view,'_input_corrected_scaled.csv'))%>%
    column_to_rownames(colnames(.)[1])
}
prot_Y_val <- read_validation('prot') %>% get_il1b()
gex_Y_val <- read_validation('gex') %>% get_il1b()
cyt_Y_val <- read_validation('cyt') %>% t()

Y_val <- cbind(data.frame('IL1B_protein' = unname(prot_Y_val), 'IL1B_gene' = unname(gex_Y_val)),cyt_Y_val)

#Predictors
covar_X_val <- fread(paste0(outdir,'validation/validation_covs.csv'))%>%
  dplyr::select(sample,AGE,SEX_BIRTH,BMI_BASELINE,season_sin,season_cos) %>%
  column_to_rownames('sample')
factor_X_val <- fread(paste0(outdir,'validation/validation_factors_alpha0.01_corrected_scaled.csv'))%>%
  column_to_rownames(colnames(.)[1])

X_val <- cbind(covar_X_val,factor_X_val)

#Modelling
explain_variance <- function(mat1, mat2){
  iterate_X <- lapply(seq_along(colnames(mat1)), function(i) colnames(mat1)[1:i])
  models <- lapply(colnames(mat2), function(y){
    
    lapply(iterate_X, function(cols){
      df_tmp <- cbind(mat1[,cols], data.frame('yc' = mat2[,y]))
      r2 <- summary(lm(yc ~ ., df_tmp ))$adj.r.squared
      return(data.frame('r2' = r2, 'predictors' = paste(cols, collapse = ' + '), 'lastpred' = cols[length(cols)]))
    })%>% bind_rows()%>%
      mutate('addr2' = r2-lag(r2, default = 0), 'addr2' = ifelse(addr2<0,0,addr2), 'y' = y)
    
  }) %>% bind_rows() %>% mutate(lastpred = factor(lastpred, levels = rev(unique(lastpred))))%>%
    group_by(y) %>% mutate('maxr2' = max(r2))
}

models_discovery <- explain_variance(X,apply(Y,2,scale))
models_validation <- explain_variance(X_val,apply(Y_val,2,scale))

models <- rbind(models_discovery %>% mutate('cohort' = 'discovery'),
                models_validation %>% mutate('cohort' = 'validation'))
fwrite(models, paste0(outdir,'figure2/explained_variance_all.csv'))

#genetics
#all_gen <- fread(cmd = 'grep -i il1b 
#                  meta_cQTL/out/2000HIV-EU-discovery/*/mapping/main_genomewide_loci.tsv', fill=TRUE)


pred_colors <- c('AGE' = 'black', 'SEX_BIRTH' = 'black', 'BMI_BASELINE' = 'black',
                 'season_sin' = 'grey', 'season_cos' = 'grey')

nfact <- ncol(get_factors(model)[[1]])
factor_colors <- rep('#FFCCCC', nfact)
names(factor_colors) <- paste0('Factor',1:nfact)
factor_colors['Factor6'] <- '#800080'

#levels
levels_y <- models %>% filter(cohort == 'validation') %>% arrange(maxr2) %>% pull(y) %>% unique()

models %>% 
  mutate(y=factor(y,levels=rev(levels_y)), cohort = factor(cohort, levels = c('discovery','validation')))%>%
  ggplot()+
  geom_col(aes(x=cohort,y=addr2,fill=lastpred,color=lastpred), position = 'stack',width = .7)+
  coord_flip()+
  facet_grid(y~.,switch = 'y')+
  scale_x_discrete(position = 'top')+
  scale_color_manual(values = c(pred_colors,factor_colors))+
  scale_fill_manual(values = c(pred_colors,factor_colors))+
  theme_bw()+theme(panel.grid = element_blank(),legend.position = 0,panel.border = element_blank(),
                   axis.line.x = element_line(),strip.text.y.left = element_text(angle = 0, hjust=1),
                   strip.background = element_blank(),axis.ticks.y = element_blank(),
                   panel.spacing = unit(0,'lines'),plot.margin = unit(c(0,0,0,0),'cm'))+
  xlab('')+
  ylab('Explained variance')

ggsave(paste0(outdir,'figure2/explained_variance_all',append,'.pdf'), width = 6, height = 14)

models %>%filter(grepl('il1b', y, ignore.case = T))%>%
  mutate(y=factor(y,levels=rev(levels_y)), cohort = factor(cohort, levels = c('discovery','validation')))%>%
  ggplot()+
  geom_col(aes(x=cohort,y=addr2,fill=lastpred,color=lastpred), position = 'stack',width = .7)+
  coord_flip()+
  facet_grid(y~.,switch = 'y')+
  scale_x_discrete(position = 'top')+
  scale_color_manual(values = c(pred_colors,factor_colors))+
  scale_fill_manual(values = c(pred_colors,factor_colors))+
  theme_bw()+theme(panel.grid = element_blank(),panel.border = element_blank(),
                   axis.line.x = element_line(),strip.text.y.left = element_text(angle = 0, hjust=1),
                   strip.background = element_blank(),axis.ticks.y = element_blank(),
                   panel.spacing = unit(0,'lines'),plot.margin = unit(c(0,0,0,0),'cm'))+
  xlab('')+
  ylab('Explained variance')
ggsave(paste0(outdir,'figure2/explained_variance_il1b',append,'.pdf'), width = 10, height = 4)

ggboxplot(distinct(models[,c('cohort','y', 'maxr2')]), x = 'cohort', y = 'maxr2', add = 'jitter',
          add.params = list(alpha=.05))+
  stat_compare_means()

ggpaired(distinct(models[,c('cohort','y', 'lastpred', 'addr2')]), x = 'cohort', y = 'addr2', id = 'lastpred', add = 'jitter',
          add.params = list(alpha=.05))+
  facet_wrap(~y)+
  stat_compare_means(paired = T,label.y = 0.2)

models%>%
  #filter(grepl('il1b', y, ignore.case = T))%>%
  .[,c('cohort','y', 'lastpred', 'addr2')]%>%
  distinct() %>%
  pivot_wider(names_from = 'cohort', values_from = 'addr2')%>%
  mutate('changer2' = validation-discovery)%>%
  group_by(lastpred)%>%mutate(medchange = median(changer2))%>%ungroup()%>%
  arrange(medchange)%>%mutate(lastpred = factor(lastpred, levels = rev(unique(lastpred))))%>%
  ggboxplot(x='lastpred',y='changer2')+
  geom_hline(aes(yintercept = 0))+
  theme(axis.text.x = element_text(angle=90, hjust = 1))

