suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
})

outdir = "/vol/projects/CIIM/2000HIV/cQTL/mofa/out/"
args = commandArgs(trailingOnly=TRUE)
corrected = args[1]
model <- readRDS(model, paste0(outdir,'model', corrected, '.rds'))

#Which factors or covariates differentiate the clusters?
clusterTest <- function(df, cluster, variable, g1=NULL, g2=NULL){
  df_tmp <- df[,c(cluster,variable)] %>% na.omit()
  colnames(df_tmp) <- c('clust','var')

  #Subset to the groups we want to test
  if(!is.null(g1) && !is.null(g2)){
    df_tmp <- df_tmp %>% filter(clust %in% c(g1,g2))
  }else if(!is.null(g1)){
    g2 <- 'NULL'
    df_tmp <- df_tmp %>% mutate(clust = ifelse(clust != g1, 'rest', clust))
  }else{
    g1 <- 'NULL'
    g2 <- 'NULL'
  }

  dif <- NA
  if(length(unique(df_tmp$var)) < 10){
    #For categorical variables, just check if they are associated
    #Quick fix, if there are less than 10 unique values I assume its categorical
    t <- table(df_tmp$clust, df_tmp$var)
    #Columns with less than 5 counts
    low_col <- apply(t, 2, function(x) any(x<5))
    #Rows
    low_row <- apply(t, 1, function(x) any(x<5))
    if(sum(low_col) < sum(low_row)){
      t <- t[,!low_col]  #Counts lower than 5 are not suitable for chisq
    }else{
      t <- t[!low_row,]  #Counts lower than 5 are not suitable for chisq

    }
    p <- chisq.test(t)$p.value
    test <- 'chisq'

  }else if(length(unique(df_tmp$clust)) > 2){
    #For multiple group comparisons
    df_tmp <- mutate_at(df_tmp, 'var', as.numeric)
    p <- kruskal.test(as.formula('var ~ clust'), df_tmp)$p.value
    test <- 'kruskal'

  }else if(length(unique(df_tmp$clust)) == 2){
    #For binary group comparisons
    df_tmp <- mutate_at(df_tmp, 'var', as.numeric)
    p <- wilcox.test(as.formula('var ~ clust'), df_tmp)$p.value
    dif <- df_tmp %>% group_by(clust) %>% summarise(var = median(var)) %>% pull(var, name = clust)
    if(g1 == 'NULL'){
      dif <- dif[2]-dif[1]
    }else if(g2 == 'NULL'){
      dif <- dif[as.character(g1)]-dif['rest']
    }else{
      dif <- dif[as.character(g1)]-dif[as.character(g2)]
    }
    test <- 'wilcox'

  }
  return(data.frame('variable' = variable, 'cluster' = cluster , 'test' = test, 'dif' = unname(dif), 'p' = p, 'g1' = g1, 'g2' = g2))
}

model_cluster_test <- lapply(colnames(model_df)[-c(1,2)] %>%
                               .[!grepl('^cluster|^UMAP|^ETHNICITY$|^EC_combined$|^ethnicity$', .)], function(c){
  print(c)
  lapply(1:4, function(g) clusterTest(model_df, cluster = 'cluster4', variable = c, g1 = g))
})%>%bind_rows()%>%
  mutate(variable_type = ifelse(grepl('^single_group.Factor', variable), 'Factor','Covariate'))%>%
  group_by(variable_type)%>%mutate('padj' = p.adjust(p))%>%
  arrange(padj)

model_df_views <- get_data(model)%>%lapply(function(x) as.data.frame(x[[1]]))%>%bind_rows()%>%
  t() %>% as.data.frame() %>% rownames_to_column('sample') %>%
  left_join(model_df %>% select(sample,cluster4))

model_cluster_test_views <- lapply(colnames(model_df_views)[-1]%>% .[!grepl('^cluster|^UMAP', .)], function(c){
  print(c)
  lapply(1:4, function(g) clusterTest(model_df_views, cluster = 'cluster4', variable = c, g1 = g))
})%>%bind_rows()%>%mutate('padj' = p.adjust(p)) %>% arrange(padj)%>%
  mutate('view' = ifelse(variable %in% rownames(model@data$cyt[[1]]), 'cyt',
                         ifelse(variable %in% rownames(model@data$metab[[1]]), 'metab', 'prot')))

ggplot(model_cluster_test_views, aes(x=view, color=view, y=sign(dif)*-log10(padj)))+
  geom_boxplot()+facet_wrap(~g1)+
  theme_bw()+theme(panel.grid = element_blank())
ggsave(paste0(outdir,'DE_views',corrected,'.pdf'))


ggplot(model_cluster_test_views %>% filter(view=='prot')%>%
         mutate(variable = gsub('_II','',variable))%>%
         rowwise()%>%
         mutate(prot_group = str_split(variable, pattern = '_', simplify = T) %>% .[length(.)]),

       aes(x=prot_group, y=sign(dif)*-log10(padj)))+
  geom_boxplot()+facet_wrap(~g1)+
  theme_bw()+theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90),legend.position = 0)
ggsave(paste0(outdir,'DE_protein_categories',corrected,'.pdf'))


ggplot(model_cluster_test_views %>% filter(view == 'cyt')%>%
         separate(variable, into = c('pbmc','time')), aes(x=time,  y=dif))+
  geom_boxplot()+
  facet_wrap(~g1)+
  geom_hline(aes(yintercept=0),linetype='dashed')+
  theme_bw()+theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90))
ggsave(paste0(outdir,'DE_cytokine_times',corrected,'.pdf'))


plot_dimred(model, method = 'UMAP', color_by = model_cluster_test_views %>% filter(view == 'prot') %>% pull(variable) %>% unique() %>% .[1:12])
ggsave(paste0(outdir,'umap_DE_prot',corrected,'.pdf'))
plot_dimred(model, method = 'UMAP', color_by = model_cluster_test_views %>% filter(view == 'cyt') %>% pull(variable) %>% unique() %>% .[1:12])
ggsave(paste0(outdir,'umap_DE_cyt',corrected,'.pdf'))
plot_dimred(model, method = 'UMAP', color_by = model_cluster_test_views %>% filter(view == 'metab') %>% pull(variable) %>% unique() %>% .[1:12])
ggsave(paste0(outdir,'umap_DE_metab',corrected,'.pdf'))
