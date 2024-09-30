suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(MOFAdata)
  library('org.Hs.eg.db')
  library(AnnotationDbi)
  library(ComplexHeatmap)
  library(ggpubr)
  library(clusterProfiler)
  library("AnnotationDbi")
  library(rstatix)
})

covs = '/vol/projects/CIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "/vol/projects/CIIM/2000HIV/cQTL/mofa/out/"

append = '_corrected_scaled'
dir.create(paste0(outdir,'profiles'))
#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Colors
source('/vol/projects/CIIM/2000HIV/cQTL/mofa/code/peRReo.R')
viewcol <- latin_palette('buenavista', n=5)
names(viewcol) <- views_names(model)


#Clin
clin_df <- fread(paste0(outdir,'clin_cor/clin_df_corrected_scaled.csv')) %>% dplyr::rename('sample' = V1)

#Add factor, covs, clin and umap
model_df <- inner_join(
  samples_metadata(model),
  get_factors(model) %>% as.data.frame() %>% rownames_to_column('sample'),by='sample',
) %>% inner_join(model@dim_red$UMAP,by='sample')  %>%
  inner_join(clin_df,by='sample')

#Taking care of the duplicated inr
all(na.omit(model_df$INR_YESNO.x == model_df$INR_YESNO.y))

model_df <- model_df %>% dplyr::rename(INR_YESNO = INR_YESNO.x ) %>%
  dplyr::select(-INR_YESNO.y)

#Plotting the umap
umap_covar <- function(df,covar,col_fun, save = T){
  p <- ggplot(df)+
    geom_point(aes(x=UMAP1,y=UMAP2,fill=.data[[covar]]),shape=21,size=2)+
    theme_bw()+xlab('')+ylab('')+
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(),panel.border = element_blank())+col_fun
  covar <- gsub('`', '', covar) %>% gsub('%', '_', .)
  legend <- cowplot::get_legend(p)
  if(save) {
    ggsave(paste0(outdir,'umaps/umap_',covar,append,'.pdf'),p+theme(legend.position = 0),width=6,height=6)
    ggsave(paste0(outdir,'umaps/umap_',covar,append,'_legend.pdf'),as_ggplot(legend),width=2,height=2)
  }
  return(p)
}

#Function to plot a feature per group
plot_data_boxplot <- function(mofa, df, feature, covar, save = F, name = ''){
  n <- sum(lengths(feature))
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n/ncol)

  grou_length <- length(unique(df[,covar]))
  newcovar <- gsub('`', '', covar) %>% gsub('%', '_', .)
  df[,newcovar] <- df[,covar]
  f <- feature
  data <- get_data(mofa,features = f, as.data.frame = T)
  pheno <- df
  df <- left_join(data,pheno,by='sample')

  df <- filter(df, !is.na(df[,newcovar]))

  if(is.numeric(df[,newcovar]) & grou_length > 10){
    p <-  ggplot(df, aes(x=.data[[newcovar]],y=value))+
      geom_point(alpha=.4)+geom_smooth(method = 'lm')+
      facet_wrap(c('feature'),ncol = ncol)+
      ylab('Scaled value')+
      stat_cor(method = 'spearman')+theme_bw()+theme(panel.grid = element_blank(),strip.background = element_blank())
  }else{
    p <- ggboxplot(df, x=newcovar,y='value', add = 'jitter', add.params = list(alpha=.1,size=.5))+
      stat_pwc(label = 'p.signif', hide.ns = T)+
      facet_wrap(c('feature'),ncol = ncol)+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
      ylab('Scaled value')+
      theme(axis.text.x = element_text(angle=90),strip.background = element_blank(),text = element_text(size=15))
  }


  if(save)ggsave(paste0(outdir,'profiles/',names(feature),'_',newcovar, name, '_dp',append,'.pdf'),
                 p,width=min(ncol,2),height=nrow*2)

  return(p)
}

#Testing the significance per group
test_signif <- function(mofa, df, feature, covar){
  grou_length <- length(unique(df[,covar]))
  newcovar <- gsub('`', '', covar) %>% gsub('%', '_', .)
  df[,newcovar] <- df[,covar]
  f <- feature
  data <- get_data(mofa,features = f, as.data.frame = T)
  pheno <- df
  df <- left_join(data,pheno,by='sample')

  df <- filter(df, !is.na(df[,newcovar]))

  if(is.numeric(df[,newcovar]) & grou_length > 10){
    test <- group_by(df,feature) %>% cor_test(vars = newcovar, vars2 = 'value')%>%pull(p,name = feature)
  }else{
    test <- group_by(df,feature) %>% wilcox_test(as.formula(paste0('value ~ ',newcovar)))%>%pull(p,name = feature)
  }
  return(test)
}

#Profile
factor_profile <- function(mofa, factor, save = T, df, covar){

  #Changing protein names too
  features_names(mofa)$prot <- features_names(mofa)$prot %>% gsub('_Inflam','-Inflam', .) %>% gsub('_Cardio','-Cardio',.)%>%
    gsub('_Onco','-Onco',.)%>%gsub('_Neuro','-Neuro',.)%>%gsub('_II$','-II',.) %>%
    str_split(pattern = '-', simplify = T) %>% .[,1]

  #Changing methylation names too
  features_names(mofa)$meth  <- features_names(mofa)$meth %>% str_split(pattern = '_', simplify = T)%>% .[,2]
  features_names(mofa)$meth <- ifelse(features_names(mofa)$meth == '', 'NA', features_names(mofa)$meth)

  #Get significant weights
  weights <- get_weights(mofa, factors = factor, scale = F)

  #Z scaling for the weights
  weights <- lapply(names(weights), function(x){
    rowna <- rownames(weights[[x]])
    df <- weights[[x]]%>%apply(2, scale)
    df <- as.data.frame(df)
    df$view <- x
    df$feature <- rowna
    colnames(df) <- c('factorx','view','feature')
    return(df)
  }
  )%>%bind_rows()

  #The sig weights
  alpha <- 0.01
  sigweights <- weights %>% filter(abs(factorx) > qnorm(alpha/2, lower.tail = F))
  fwrite(sigweights,paste0(outdir,'profiles/Sigweights_',factor,append,'.csv'))

  #Plot the sigweights: Top 10 per view
  sigweights <- group_by(sigweights, view)

  lapply(unique(sigweights$view), function(x) {
    y <- sigweights %>% slice_max(abs(factorx), n = 10) %>% filter(view == x) %>% pull(feature, name = view)
    y <- split(y, names(y))
    print(y)
    plot_data_boxplot(mofa = mofa, df = df, feature = y,covar = covar, save = T, name = paste0('_factor', factor))
    plot_data_boxplot(mofa = mofa, df = df, feature = y,covar = paste0('single_group.Factor',factor), save = T)
    sig_y <- sigweights %>% filter(feature!='NA')%>% slice_max(abs(factorx), n = 50) %>% filter(view == x) %>% pull(feature, name = view)
    sig_y <- split(sig_y, names(sig_y))
    sig_y <- test_signif(mofa = mofa,df = df,feature = sig_y,covar = covar)
    if(any(sig_y < 0.05)){
      sig_y <- sig_y %>% .[.<0.05] %>% .[1:min(2,length(.))]
      sig_y <- list(names(sig_y))
      names(sig_y) <- names(y)
      plot_data_boxplot(mofa = mofa, df = df, feature = sig_y,covar = covar,
                        save = T, name = paste0('_factor', factor,'_significant'))
    }
  })

  #Plot those that are significant for the covar too


  #Any weights with common sign between gex and prot?
  weights_gex_prot <- weights %>% filter(view %in% c('prot','gex')) %>%
    group_by(feature)%>%mutate(n=n_distinct(view))%>%filter(n==2)%>%
    pivot_wider(names_from = 'view', values_from = 'factorx',values_fn = function(x)x[1])%>%
    mutate(combined=gex*prot) %>%
    filter(combined > qnorm(alpha, lower.tail = F)*qnorm(alpha, lower.tail = F))

  if(length(weights_gex_prot$feature) > 0){
    plot_data_boxplot(mofa = mofa, df = df, feature = list('gex'=weights_gex_prot$feature),covar = covar, save = T, name = paste0('_factor', factor,'_gex_prot'))
    plot_data_boxplot(mofa = mofa, df = df, feature = list('prot'=weights_gex_prot$feature),covar = covar, save = T, name = paste0('_factor', factor,'_gex_prot'))
    plot_data_boxplot(mofa = mofa, df = df, feature = list('gex'=weights_gex_prot$feature),covar = paste0('single_group.Factor',factor), save = T,name = '_gex_prot')
    plot_data_boxplot(mofa = mofa, df = df, feature = list('prot'=weights_gex_prot$feature),covar = paste0('single_group.Factor',factor), save = T,name = '_gex_prot')

  }



  #Enrichments
  btm <- read.gmt(paste0(outdir,'BTM_for_GSEA_20131008.gmt')) %>%
    #mutate(gene2 = bitr(gene,fromType = 'SYMBOL', toType = 'ENSEMBL', drop = F, OrgDb="org.Hs.eg.db"))
    mutate(gene2 = 1) %>% filter(!grepl('^TBA ', term))%>%
    pivot_wider(names_from = gene, values_from = gene2, values_fill = 0) %>%
    column_to_rownames('term') %>% as.matrix()

  #My own cytokine gmt
  cytgmt <-  features_names(model)$cyt %>% str_split(pattern = '_', simplify = T) %>% as.data.frame() %>%
    mutate('a' = paste(V1,V2,V3,V4, sep = '_')) %>% pivot_longer(col = c(V1,V2,V3,V4))%>%
    dplyr::select(-name)%>%mutate('c' = 1)%>%filter(value != 'pbmc')%>%
    pivot_wider(names_from = a, values_from = c, values_fill = 0) %>%
    column_to_rownames('value')%>% as.matrix()

  enr_view <- function(mofa, factor, view, gmt, bysign = T){
    if(bysign){
      enr_pos <- run_enrichment(mofa, view = view, feature.sets = gmt, factors = factor, sign = 'pos', min.size = 5)
      enr_neg <- run_enrichment(mofa, view = view, feature.sets = gmt, factors = factor, sign = 'neg', min.size = 5)

      enr <- rbind(
        as.data.frame(enr_pos$set.statistics) %>% rownames_to_column('path') %>%
          filter(path %in% enr_pos$sigPathways[[1]]) %>% mutate('sign' = 1),
        as.data.frame(enr_neg$set.statistics) %>% rownames_to_column('path') %>%
          filter(path %in% enr_neg$sigPathways[[1]]) %>% mutate('sign' = -1)
      )%>%dplyr::rename('statistic' = colnames(.)[2]) %>% mutate(view = view)
    }else{
      enr_all <- run_enrichment(mofa, view = view, feature.sets = gmt, factors = factor, min.size = 5)
      enr <-  as.data.frame(enr_all$set.statistics) %>% rownames_to_column('path') %>%
        filter(path %in% enr_all$sigPathways[[1]]) %>% mutate('sign' = 1)%>%
        dplyr::rename('statistic' = colnames(.)[2]) %>% mutate(view = view)
    }
    #if(nrow(enr) == 0){
    #  enr <- data.frame('path' = NA, 'statistic' = NA, 'sign' = NA, 'view' = view)
    #}


  }

  #Enrichment for prot, gex and cyt
  enr_all <- rbind(
    enr_view(mofa, factor, 'prot', btm),
    enr_view(mofa, factor, 'gex', btm),
    enr_view(mofa, factor, 'cyt', cytgmt) #,
    #enr_view(mofa, factor, 'meth', btm, bysign = T)
  )
  fwrite(enr_all,paste0(outdir,'profiles/Enrichment_',factor,append,'.csv'))

  enr_all <- enr_all %>%
    group_by(view) %>% slice_max(abs(statistic),n=20)%>%
    mutate(path = factor(path, levels =path,labels = str_split(path,'\\(M', simplify = T)[,1] ))
                           #gsub('(.{1,30})(\\s|$)', '\\1\n', .)))

  changeIds <- function(ids){
    mapIds(x = org.Hs.eg.db,
           keys = ids,
           column = "SYMBOL",
           keytype = "ENSEMBL",
           multiVals = "first") %>% as.data.frame() %>% rownames_to_column('ENSEMBL')%>%
      dplyr::rename('SYMBOL' = colnames(.)[2])%>%
      mutate(SYMBOL = ifelse(is.na(SYMBOL), ENSEMBL,SYMBOL))%>%pull(SYMBOL)
  }

  data("reactomeGS")
  colnames(reactomeGS) <- changeIds( colnames(reactomeGS))
  data("MSigDB_v6.0_C2_human")
  colnames(MSigDB_v6.0_C2_human) <- changeIds(colnames(MSigDB_v6.0_C2_human))
  data("MSigDB_v6.0_C5_human")
  colnames(MSigDB_v6.0_C5_human) <- changeIds(colnames(MSigDB_v6.0_C5_human))

  enr_reactome <- rbind(enr_view(mofa,factor,'prot',reactomeGS),enr_view(mofa,factor,'gex',reactomeGS))%>%
    slice_max(abs(statistic),n=20)
  #enr_msig <- rbind(enr_view(mofa,factor,'prot',MSigDB_v6.0_C2_human),enr_view(mofa,factor,'gex',MSigDB_v6.0_C2_human))
  enr_go <- rbind(enr_view(mofa,factor,'prot',MSigDB_v6.0_C5_human),enr_view(mofa,factor,'gex',MSigDB_v6.0_C5_human))%>%
    slice_max(abs(statistic),n=20)


  #Metabolite enrichment (overrepresentation)
  sigmetab <- sigweights %>% filter(view == 'metab') %>% pull(feature)
  metab_ids  <- fread('/vol/projects/CIIM/2000HIV/Metabolites/metabo_index.tsv') %>%
    filter(name %in% sigmetab)
  #Jianbo's approach
  sigmetab_ids <- metab_ids$HMDB_ids %>% str_split(',', simplify = T) %>% .[,1]
  #Xun's approach
  #sigmetab_ids <- metab_ids$HMDB_ids %>% str_split(',', simplify = T) %>% as.vector() %>% .[.!='']

  #Barplot
  Molecular_profile_plot <- function(enr_df){
    enr_df %>%
      ggplot(aes(x=reorder(path, sign*statistic), y=sign*statistic, fill=view))+
      geom_col()+
      scale_fill_manual(values = viewcol)+
      #geom_text(data = enr_all %>% filter(sign==1), aes(label=path), angle = 90, size = 3, hjust=1,vjust=.5)+
      #geom_text(data = enr_all %>% filter(sign==-1), aes(label=path), angle = 90, size = 3, hjust=0,vjust=.5)+
      facet_grid(~view, scales = 'free_x', space = 'free')+
      theme_bw()+theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90,hjust=1,vjust=.75),
                       legend.position = 0,axis.ticks.x = element_blank(),)+
      theme(strip.background = element_blank(),
            plot.margin = unit(c(0,0,0,0), 'cm'),strip.placement = 'outside')+
      xlab('')+ylab(paste('Factor',factor))
  }

  mini_Molecular_profile_plot <- function(enr_df){
    enr_df_pos <- enr_df %>% filter(sign ==1) 
    enr_df_neg <- enr_df %>% filter(sign == -1) 
    if(!nrow(enr_df_pos)>1){
      enr_df <- enr_df_neg
    }else if(!nrow(enr_df_neg)>1){
      enr_df <- enr_df_pos
    }else{
      enr_df <- rbind(enr_df_pos %>% group_by(view) %>% slice_max(statistic),
                      enr_df_neg %>% group_by(view) %>% slice_max(statistic))
    }
    
    enr_df %>%
      ggplot(aes(x=reorder(path, sign*statistic), y=sign, fill=view))+
      geom_col()+
      geom_text(data = enr_df,aes(y=sign/2,label=path), angle = 90, size = 2, hjust=0.5,vjust=.5)+
      #geom_text(data = enr_df %>% filter(sign==-1),
      #aes(y=sign,label=path), angle = 90, size = 2,hjust=0,vjust=.5)+
      scale_fill_manual(values = viewcol)+
      facet_grid(~view, scales = 'free_x', space = 'free')+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0)))+
      theme_bw()+theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90,hjust=1,vjust=.75),
                       legend.position = 0,axis.ticks = element_blank(),axis.text.y = element_blank(),
                       axis.line.y.left = element_line(), panel.border = element_blank())+
      theme(strip.background = element_blank(),axis.text.x = element_blank(),
            plot.margin = unit(c(0,0,0,0), 'cm'),strip.placement = 'outside')+
      xlab('')+ylab(paste('Factor',factor))
  }


  Molecular_profile_enrichment <- Molecular_profile_plot(enr_all)
  mini_Molecular_profile_enrichment <- mini_Molecular_profile_plot(enr_all)
  Molecular_profile_reactome <- Molecular_profile_plot(enr_reactome)
  Molecular_profile_go <- Molecular_profile_plot(enr_go)


  if(save){
    ggsave(paste0(outdir,'profiles/Molecular_profile_enrichment_factor',factor,append,'.pdf'),
                 Molecular_profile_enrichment,width=ceiling(nrow(enr_all)/8),height=6)
    ggsave(paste0(outdir,'profiles/mini_Molecular_profile_enrichment_factor',factor,append,'.pdf'),
                  mini_Molecular_profile_enrichment,width=length(unique(enr_all$view)),height=3)
    ggsave(paste0(outdir,'profiles/Molecular_profile_reactome_factor',factor,append,'.pdf'),
           Molecular_profile_reactome,width=ceiling(nrow(enr_reactome)/8),height=8)
    ggsave(paste0(outdir,'profiles/Molecular_profile_go_factor',factor,append,'.pdf'),
           Molecular_profile_go,width=ceiling(nrow(enr_go)/8),height=8)
  }


}



#Plotting UMAP of the covar, UMAP of the factor, boxplot or correlation
full_factor_profile <- function(mofa, df, covar, factor, col.fun = NULL, col.values = NULL){

  #Prepare things
  full_factor <- paste0('single_group.Factor',factor)
  grou_length <- length(unique(df[,covar]))

  #Fix weird covars
  df[,gsub('%','_', covar)] <- df[,covar]
  covar <- gsub('%','_', covar)


  #UMAP covar
  if(!is.null(col.fun)){
    umap_covar(df, covar, col_fun = col.fun)
  }else if(!is.null(col.values)){
    umap_covar(df, covar, scale_fill_manual(values = col.values))
  }else if(is.numeric(df[,covar]) & grou_length > 10){
    umap_covar(df, covar, scale_fill_gradient2(low = 'darkblue', mid = 'white',high = 'darkred', midpoint = median(df[,covar])))
  }else{
    umap_covar(df %>% mutate_at(covar, as.character), covar, scale_fill_discrete())
  }

  #UMAP factor
  umap_covar(df, full_factor, scale_fill_gradient2(low = 'darkblue', mid = 'white',high = 'darkred'))

  #Confirm factor to covar
  if(is.numeric(df[,covar]) & grou_length > 10){
    p <-  ggplot(df, aes(x=.data[[covar]],y=.data[[full_factor]]))+
      geom_point(alpha=.4)+geom_smooth(method = 'lm')+
      stat_cor()+theme_bw()+theme(panel.grid = element_blank())+ylab(paste('Factor', factor))
  }else{
    p <- ggboxplot(filter(df, !is.na(df[,covar])), x=covar,y=full_factor)+
    ylab(paste('Factor', factor))+
      stat_pwc(label = 'p.signif', hide.ns = T)+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
      theme(axis.text.x = element_text(angle=90),text = element_text(size=15))
  }
  ggsave(paste0(outdir,'profiles/Factor_',factor,'_',covar,'_dp',append,'.pdf'), plot = p,width=2,height=2)

  #factor_profile
  factor_profile(mofa,factor, df = df, covar = covar, save = T)
}

full_factor_profile(model, model_df,'PLAQUE_YESNO',6)
full_factor_profile(model, model_df,'days_post_lockdown',1)
full_factor_profile(model, model_df,'MH_TR_CV',8)
full_factor_profile(model, model_df,'MH_TR_CV_D%Myocardial_infarction',8)
full_factor_profile(model, model_df,'MH_TR_RESP_D%COPD',11)
full_factor_profile(model, model_df,'MH_TR_RESP_D%COPD',11)
full_factor_profile(model, model_df,'RP_YESNO',20)


#Adding EEAA
eaa <- fread('/vol/projects/CIIM/2000HIV/AnalysisDV/Discovery_MethylationDataBased/3Aging/Discovery.DNAmAge.output.csv')%>%
  dplyr::select(c("SampleID","EEAA"))%>%
  dplyr::rename(sample=SampleID)

model_df <- left_join(model_df,eaa, by = 'sample')
full_factor_profile(model, model_df,'EEAA', 20)
full_factor_profile(model, model_df,'EEAA', 4)



#All profiles
lapply(1:length(colnames(clin_df)[grepl('Factor',colnames(clin_df))]),
       function(f) {
         print(f)
         capture.output(factor_profile(model, f, T, model_df, 'MH_TR_CV_D%Myocardial_infarction'))}
)

