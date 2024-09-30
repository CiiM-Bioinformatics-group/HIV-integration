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
dir.create(paste0(outdir,'umaps'))
dir.create(paste0(outdir,'overview'))

for(i in c('corrected','uncorrected')){
  for(j in c('scaled','unscaled')){
    append = paste0('_',i,'_',j)
    print(append)

    modelfile = paste0(outdir,"model",append,".hdf5")

    #Colors: From jbgb13/peRReo
    source('peRReo.R')
    viewcol <- latin_palette('buenavista', n=5)

    #Load
    model <- load_model(modelfile)

    #Remove highly correlated factors with the number of features
    factors <- do.call("rbind",get_factors(model))
    r <- suppressWarnings( t(do.call('rbind', lapply(model@data, function(x)
      abs(cor(colMeans(do.call("cbind",x),na.rm=TRUE),factors, use="pairwise.complete.obs"))
    ))) )
    colnames(r) <- views_names(model)

    pdf(paste0(outdir, 'overview/correlation_to_average',append,'.pdf'),width = 4,height = 8)
    Heatmap(r,cluster_rows = F,cluster_columns = F)
    dev.off()

    #Very strict, only those lower than 0.6
    intercept_factors <- which(apply(r,1,function(x)any(x>0.6)))
    goodfactors <- rownames(r)[-intercept_factors]

    #First of all, noncorrelated factors, so variance is recalculated
    model <- subset_factors(model, goodfactors)

    #topfactors
    topfactors <- model@cache$variance_explained$r2_per_factor[[1]] %>%
      apply(1,sum)
    topfactors <- names(topfactors)[topfactors > 1]

    #Only the good and the top
    model <- subset_factors(model, topfactors)

    #Changing names

    #the gene names
    #New gene names (Maybe add this to mofa_load)
    genes_key <- features_names(model)$gex %>%
      str_split(pattern = '\\.', simplify = T) %>% .[,1] %>%
      bitr(.,toType = 'SYMBOL', fromType = 'ENSEMBL', drop = F, OrgDb="org.Hs.eg.db") %>%
      pull(SYMBOL, name = ENSEMBL)
    genes_key[is.na(genes_key)] <- names(genes_key[is.na(genes_key)])

    features_names(model)$gex  <- genes_key[
      features_names(model)$gex %>%
        str_split(pattern = '\\.', simplify = T) %>% .[,1] ]

    #Metabolite names
    metab_ids  <- fread('2000HIV/Metabolites/ionMz2ID.tsv')%>%
                  pull(name, name=formula)
    features_names(model)$metab <- metab_ids[features_names(model)$metab]

    #Methylation names
    meth_df  <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)%>%
      as.data.frame()%>%dplyr::select(Name, UCSC_RefGene_Name, chr, pos)%>%
      filter(Name %in%  features_names(model)$meth)

    #How are the selected probes in the genome?
    ggplot(meth_df %>% mutate(chr = factor(chr, levels = paste0('chr',1:22))),
           aes(x=pos,y=1))+geom_jitter()+facet_grid(~chr, scales = 'free_x', space = 'free')+
      theme_bw()+
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90, hjust=1,vjust=0.2),
            panel.spacing = unit(0, "lines"), strip.background = element_blank(),
            axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(),
            axis.title.x = element_blank(), panel.border = element_blank(),axis.line.y = element_line(),
            plot.margin = unit(c(0,0,0,0), 'cm'),strip.placement = 'outside')
    ggsave(paste0(outdir, 'selected_meth_probes.pdf'))

    meth_df <- meth_df %>% separate(UCSC_RefGene_Name, into = 'UCSC_RefGene_Name', sep = ';')
    meth_ids <- meth_df %>% mutate('name_gene' = paste(Name, UCSC_RefGene_Name, sep = '_'))%>%
      pull(name_gene, name=Name)

    features_names(model)$meth <- meth_ids[features_names(model)$meth]

    #Plot
    plot_data_overview(model, colors = viewcol)
    ggsave(paste0(outdir, 'overview/data_overview',append,'.pdf'))

    #Covariates
    #Reading the matrixqtl files
    read_tsv <- function(file, name){
      table <- fread(file) %>% column_to_rownames(colnames(.)[1] )%>% t() %>% as.data.frame()
      table <- table %>% rownames_to_column('sample') %>% mutate(sample = gsub('OLVG','OLV',sample))
      table <- melt(table, variable.name='feature')
      #table$view <- name
      return(table)
    }
    covdf <- fread(covs) %>% dplyr::select(-c(DNAm_SamplePlate,WGS,Cohort,ETHNICITY.Pacific_Islander,ETHNICITY.Native_American))%>%
      dplyr::rename('sample' = Record.Id)
    covdf <- left_join(samples_metadata(model), covdf, 'sample')
    samples_metadata(model) <- covdf

    #Plots
    variance_perFactorMelted = reshape2::melt(model@cache$variance_explained$r2_per_factor)
    ggplot(variance_perFactorMelted, aes(x=Var1, y=value,fill=Var2))+geom_col(position = 'stack')+theme_bw()+
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle=90))+
      scale_fill_manual(values = viewcol)
    ggsave(paste0(outdir,'overview/variance_explained_barplot',append,'.pdf'),width = 5,height = 5)

    plot_variance_explained(model, x="view", y="factor")
    ggsave(paste0(outdir,'overview/variance_explained',append,'.pdf'),width = 5,height = 8)

    #Cluster samples
    model_clust_k <- sapply(1:15,function(x) cluster_samples(model, k = x)$tot.withinss)

    pdf(paste0(outdir,'overview/Elbow_plot',append,'.pdf'),width = 4,height = 4)
    plot(1:15,model_clust_k, type = 'b')
    dev.off()

    model_clust_2 <- cluster_samples(model,2)
    model_clust_3 <- cluster_samples(model,3)
    model_clust_4 <- cluster_samples(model,4)
    model_clust_5 <- cluster_samples(model,5)
    model_clust_6 <- cluster_samples(model,6)

    #Add to the metadata
    samples_metadata(model) <- samples_metadata(model) %>%
      left_join(
        data.frame('cluster2' = model_clust_2$cluster)%>% rownames_to_column('sample'),
        by = 'sample'
      ) %>%
      left_join(
        data.frame('cluster3' = model_clust_3$cluster)%>% rownames_to_column('sample'),
        by = 'sample'
      ) %>%
      left_join(
        data.frame('cluster4' = model_clust_4$cluster)%>% rownames_to_column('sample'),
        by = 'sample'
      ) %>%
      left_join(
        data.frame('cluster5' = model_clust_5$cluster)%>% rownames_to_column('sample'),
        by = 'sample'
      ) %>%
      left_join(
        data.frame('cluster6' = model_clust_6$cluster)%>% rownames_to_column('sample'),
        by = 'sample'
      )

    #UMAP
    set.seed(1997)
    model <- run_umap(model)
    plot_dimred(model, method = 'UMAP', color_by = c('AGE','SEX_BIRTH','Institute.Abbreviation','ETHNICITY.White'))
    ggsave(paste0(outdir,'umaps/umap_covars',append,'.pdf'))

    plot_dimred(model, method = 'UMAP', color_by = c('season_cos','season_sin','cplv','days_post_lockdown'))
    ggsave(paste0(outdir,'umaps/umap_season_lockdown',append,'.pdf'))

    plot_dimred(model, method = 'UMAP', color_by = c('cluster2','cluster3','cluster4','cluster5','cluster6'))
    ggsave(paste0(outdir,'umaps/umap_cluster',append,'.pdf'))

    plot_dimred(model, method = 'UMAP', color_by = c('CD4T','CD8T','cplv'))
    ggsave(paste0(outdir,'umaps/umap_covars2',append,'.pdf'))

    model_df <- inner_join(
      samples_metadata(model),
      get_factors(model) %>% as.data.frame() %>% rownames_to_column('sample'),
    ) %>% inner_join(model@dim_red$UMAP)

    #Plotting UMAPs
    umap1 <- ggplot(model_df)+
      geom_point(aes(x=UMAP1,y=UMAP2, fill = ''),shape=21,size=2)+
      theme_bw()+xlab('')+ylab('')+scale_fill_manual(values = 'white')+
      theme(panel.grid = element_blank(), axis.ticks = element_blank(),
            axis.text = element_blank(),panel.border = element_blank(), legend.position = 0)
    ggsave(paste0(outdir,'umaps/umap_no_color',append,'.pdf'),umap1,width=6,height=6)

    umap_covar <- function(df,covar,col_fun, save = T){
      p <- ggplot(df)+
        geom_point(aes_string(x='UMAP1',y='UMAP2',fill=covar),shape=21,size=2)+
        theme_bw()+xlab('')+ylab('')+
        theme(panel.grid = element_blank(), axis.ticks = element_blank(),
              axis.text = element_blank(),panel.border = element_blank())+col_fun
      if(save)ggsave(paste0(outdir,'umaps/umap_',covar,append,'.pdf'),p+theme(legend.position = 0),width=6,height=6)
      return(p)
    }

    umap2 <- model_df %>%
      mutate('Sex' = ifelse(SEX_BIRTH ==1, 'Female', 'Male'))%>%
      umap_covar(covar = 'Sex',col_fun = scale_fill_manual(values = c('orange','darkolivegreen3'), na.value = 'White'))


    model_df <- model_df %>%
      mutate(ethnicity = ifelse(ETHNICITY.White == 1, 'White',''))%>%
      mutate(ethnicity = ifelse(ETHNICITY.Black == 1, 'Black',ethnicity))%>%
      mutate(ethnicity = ifelse(ETHNICITY.Hispanic == 1, 'Hispanic',ethnicity))%>%
      mutate(ethnicity = ifelse(ETHNICITY_MIXED == 1, 'Mixed',ethnicity))%>%
      mutate(ethnicity = ifelse(ETHNICITY.Asian == 1, 'Asian',ethnicity))%>%
      mutate(ethnicity = as.factor(ethnicity))


    ethnicity_cols <- RColorBrewer::brewer.pal(name = 'Set3', n=6)
    names(ethnicity_cols) <- rev(levels(model_df$ethnicity))
    ethnicity_cols['other'] <- 'White'

    umap3 <- umap_covar(model_df, 'ethnicity',
      scale_fill_manual(values = ethnicity_cols, na.value = 'White'))


    UMAP_ethnicity <- function(eth){
      temp_df <- model_df %>% mutate('ethnicity' = ifelse(ethnicity %in% eth, as.character(ethnicity), 'other'))
      p <- ggplot(temp_df)+
        geom_point(aes(x=UMAP1,y=UMAP2, fill= ethnicity ),shape=21,size=2)+
        labs(fill="Sex")+scale_fill_manual(values = ethnicity_cols, na.value = 'White')+
        theme_bw()+xlab('')+ylab('')+
        theme(panel.grid = element_blank(), axis.ticks = element_blank(),
              axis.text = element_blank(),panel.border = element_blank(),legend.position = 0)
      ggsave(paste0(outdir,'umaps/umap_',eth,append,'.pdf'),p,width=6,height=6)
      return(p)
    }
    umap4 <- UMAP_ethnicity(c('Black'))
    umap5 <- UMAP_ethnicity(c('Asian'))
    umap6 <- UMAP_ethnicity(c('Hispanic'))
    umap7 <- UMAP_ethnicity(c('Mixed'))

    umap8 <- umap_covar(model_df, 'AGE',
      scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = median(model_df$AGE,na.rm = T)))

    umap9 <- umap_covar(model_df, 'BMI_BASELINE',
      scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = median(model_df$BMI_BASELINE,na.rm = T)))

    umap10 <- umap_covar(model_df, 'season_sin',
      scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = median(model_df$season_sin,na.rm = T)))

    umap11 <- umap_covar(model_df, 'season_cos',
      scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = median(model_df$season_sin,na.rm = T)))

    umap12 <- umap_covar(model_df, 'cplv', scale_fill_brewer(palette = 'Accent'))

    (umap1+umap2+umap3)/(umap4+umap5+umap6)/(umap7+umap8+umap9)/(umap10+umap11+umap12)
    ggsave(paste0(outdir,'umaps/umap_all',append,'.pdf'),width=18,height=24)

    saveRDS(model, paste0(outdir,'model', append, '.rds'))
    #Save factors
    factors <- get_factors(model)[[1]] %>% t()
    fwrite(as.data.frame(factors),paste0(outdir, 'factors', append,'.txt'), sep = '\t', row.names=T)

  }
}

#Compare scaled and unscaled models
#m1 <- load_model(paste0(outdir,"model_corrected_scaled.hdf5"))
#m2 <- load_model(paste0(outdir,"model_corrected_unscaled.hdf5"))

#f1 <- get_factors(m1)[[1]]
#colnames(f1) <- paste0(colnames(f1),'_scaled')

#f2 <- get_factors(m2)[[1]]
#colnames(f2) <- paste0(colnames(f2),'_unscaled')

#Heatmap(psych::cor2(f1, f2),cluster_rows = F, cluster_columns = F)
