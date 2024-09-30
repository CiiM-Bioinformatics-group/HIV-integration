suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(ComplexHeatmap)
  library(ggpubr)
  library(parallel)
  library(rstatix)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'cov_cor'))
#mclapply(c('corrected','uncorrected'), function(i){
#mclapply(c('scaled','unscaled'), function(j){
#append = paste0('_',i,'_',j)
append = '_corrected_scaled'
    #args = commandArgs(trailingOnly=TRUE)
    #if(is.na(args[1])) corrected = '' else corrected = args[1]
    #modelfile = paste0(outdir,"model_corrected_unscaled.hdf5")

    #Load
    model <- readRDS(paste0(outdir,'model', append, '.rds'))

    #Factors
    factors <- get_factors(model, factors = "all") %>% .$single_group
    covdf <- samples_metadata(model)

    factors_cov <- inner_join(
      covdf %>% select(-group),
      factors %>% as.data.frame() %>% rownames_to_column('sample'),
      by = 'sample'
    )%>%column_to_rownames('sample')%>%
      mutate_at("CMV_IgG_IU.mL", as.numeric)
    #  mutate_at(c('cplv','Institute.Abbreviation'), function(x) as.numeric(as.factor(x)))

    #Remove intercorrelated things
    factors_cov <- factors_cov %>%
      dplyr::select(-c(cov,vac,iso,iso2, ETHNICITY,Vaccine_Group, days_post_cov, days_post_vacc,INR_YESNO,
                       colnames(.)[grepl('^cluster',colnames(.))]))
    #Transform to dummy
    factors_cov <- model.matrix(~ . -1 ,factors_cov)%>%
      as.data.frame()

    #Single correlation
    factors_cov_cor <- psych::corr.test(
      factors_cov[,grepl('^Factor',colnames(factors_cov))],
      factors_cov[,!grepl('^Factor|^cluster',colnames(factors_cov))],
      use = 'pairwise.complete')

    fwrite(as.data.frame(factors_cov_cor[c('r','p','p.adj')]), paste0(outdir,'factors_cov_cor',append,'.csv'), row.names = T)

    pdf(paste0(outdir,'cov_cor/Factors_cov_cor',append,'.pdf'),width = 10,height = 10)
    corrplot::corrplot(factors_cov_cor$r, col = colorRampPalette(colors = c('Blue', 'white','darkred'))(100),
                       tl.col = 'black',na.label = ' ', is.corr = T, addgrid.col = NA, outline = F,mar=c(0,0,2,0),
                       p.mat = factors_cov_cor$p.adj,sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
                       insig = 'label_sig', pch.col = 'black')
    dev.off()

    #Those with at least one significant hit
    onesig <- apply(factors_cov_cor$p.adj, 2, function(x) sum(x<0.05)) %>% .[.>0] %>% na.omit()

    pdf(paste0(outdir,'cov_cor/Factors_cov_cor_reduced',append,'.pdf'),width = 10,height = 10)
    corrplot::corrplot(factors_cov_cor$r[,names(onesig)], col = colorRampPalette(colors = c('Blue', 'white','darkred'))(100),
                       tl.col = 'black',na.label = ' ', is.corr = T, addgrid.col = NA, outline = F,mar=c(0,0,2,0),
                       p.mat = factors_cov_cor$p.adj[,names(onesig)],sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
                       insig = 'label_sig', pch.col = 'black')
    dev.off()

    #Testing
    factors_cov_tests <- lapply(colnames(factors_cov)[!grepl('^Factor',colnames(factors_cov))], function(var){
      print(var)

      facts <- colnames(factors_cov)[grepl('^Factor', colnames(factors_cov))]
      obs <- length(na.omit(unique(factors_cov[,var])))
      empty_df <- data.frame('p' = NA, 'factor' = NA, 'clin_var' = var, 'test' = NA, 'estimate'=NA,'group1'=NA,'group2'=NA,'n1'=NA,'n2'=NA)

      tmp.df <- factors_cov %>% dplyr::select(all_of(c(facts, var)))
      colnames(tmp.df)[ncol(tmp.df)] <- 'variable'
      tmp.df <- tmp.df %>% pivot_longer(cols = all_of(facts), names_to = 'factor', values_to = 'value')%>%
        group_by(factor)

        if(obs < 2){
          return(empty_df)
        }else if(obs== 2){
          out <- tmp.df %>%
            wilcox_test(value ~ variable, detailed = T, ref.group = '1')%>%
            dplyr::select(factor, p,estimate, group1, group2, n1, n2,conf.low,conf.high)%>%mutate('clin_var' = var, 'test' = 'wilcox')
        }else if(obs < 10){
          out <- tmp.df %>%
            kruskal_test(value ~ variable)%>%
            dplyr::select(factor, p)%>%mutate('clin_var' = var, 'test' = 'kruskal','estimate'=NA,'group1'=NA,
            'group2'=NA,'n1'=NA,'n2'=NA,'conf.low'=NA,'conf.high'=NA)
        }else if(is.numeric(factors_cov[,var])){
          out <- tmp.df %>%
            group_by(factor)%>%
            cor_test(vars = 'variable', vars2 = 'value')%>%
            dplyr::select(factor, p,cor,conf.low,conf.high)%>%mutate('clin_var' = var, 'test' = 'pearson','group1'=NA,'group2'=NA,'n1'=NA,'n2'=NA)%>%
            dplyr::rename(estimate=cor)
        }else{
          return(empty_df)
      }
    }) %>% bind_rows()%>% mutate('padj' = p.adjust(p, method = 'fdr'))
    fwrite(factors_cov_tests, paste0(outdir,'cov_cor/factors_cov_test',append,'.csv'), row.names = T)


#  }, mc.cores=2)
#}, mc.cores=4)
