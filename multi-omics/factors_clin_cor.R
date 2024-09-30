suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(ComplexHeatmap)
  library(ggpubr)
  library(rstatix)
  library(parallel)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/clin_cor/"
#mclapply(c('corrected','uncorrected'), function(i){
#  mclapply(c('scaled','unscaled'), function(j){
#    append = paste0('_',i,'_',j)
    append = '_corrected_scaled'

    #args = commandArgs(trailingOnly=TRUE)
    #corrected = args[1]

    #Read the clean table created by clin_df.R
    clin_df_clean <- fread(paste0(outdir,'clin_df',append,'.csv'))%>%
      column_to_rownames(colnames(.)[1])

    #Just correlation
    clin_df_cor <- psych::corr.test(
      clin_df_clean[,grepl('^Factor',colnames(clin_df_clean))],
      clin_df_clean[,!grepl('^Factor|^cluster',colnames(clin_df_clean))],
      use = 'pairwise.complete')

    fwrite(as.data.frame(clin_df_cor[c('r','p','p.adj')]), paste0(outdir,'clin_df_cor',append,'.csv'), row.names = T)


    pdf(paste0(outdir,'Factors_clin_cor',append,'.pdf'),width = 22,height = 10)
    corrplot::corrplot(clin_df_cor$r, col = colorRampPalette(colors = c('Blue', 'white','darkred'))(100),
                       tl.col = 'black',na.label = ' ', is.corr = T, addgrid.col = NA, outline = F,mar=c(0,0,2,0),
                       p.mat = clin_df_cor$p.adj,sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
                       insig = 'label_sig', pch.col = 'black')
    dev.off()

    #Those with at least one significant hit
    onesig <- apply(clin_df_cor$p.adj, 2, function(x) sum(x<0.05)) %>% .[.>0] %>% na.omit()

    pdf(paste0(outdir,'Factors_clin_cor_reduced',append,'.pdf'),width = 10,height = 10)
    corrplot::corrplot(clin_df_cor$r[,names(onesig)], col = colorRampPalette(colors = c('Blue', 'white','darkred'))(100),
                       tl.col = 'black',na.label = ' ', is.corr = T, addgrid.col = NA, outline = F,mar=c(0,0,2,0),
                       p.mat = clin_df_cor$p.adj[,names(onesig)],sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
                       insig = 'label_sig', pch.col = 'black')
    dev.off()

    #Wilcoxon for those with two factors, kruskal for the rest
    clin_df_tests <- lapply(colnames(clin_df_clean)[!grepl('^Factor|^cluster',colnames(clin_df_clean))], function(var){
      print(var)

      facts <- colnames(clin_df_clean)[grepl('^Factor', colnames(clin_df_clean))]
      obs <- length(na.omit(unique(clin_df_clean[,var])))
      empty_df <- data.frame('p' = NA, 'factor' = NA, 'clin_var' = var, 'test' = NA, 'estimate'=NA,'group1'=NA,'group2'=NA,'n1'=NA,'n2'=NA)

      tmp.df <- clin_df_clean %>% dplyr::select(all_of(c(facts, var)))
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
      }else if(is.numeric(clin_df_clean[,var])){
        out <- tmp.df %>%
          cor_test(value ~ variable)%>%
          dplyr::select(factor, p,cor,conf.low,conf.high)%>%mutate('clin_var' = var, 'test' = 'pearson','group1'=NA,'group2'=NA,'n1'=NA,'n2'=NA)%>%
          dplyr::rename(estimate=cor)
      }else{
        return(empty_df)
        }
    }) %>% bind_rows()%>% mutate('padj' = p.adjust(p, method = 'fdr'))
    fwrite(clin_df_tests, paste0(outdir,'clin_df_test',append,'.csv'), row.names = T)


    #Correct
    clin_df_tests_pmat <- clin_df_tests %>%
      dplyr::select(factor, clin_var, padj)%>%
      pivot_wider(names_from = clin_var,values_from = padj)%>%
      filter(!is.na(factor))%>%
       column_to_rownames('factor')%>%as.matrix()

    clin_df_tests_emat <- clin_df_tests %>%
      mutate(estimate = ifelse(is.na(estimate),1,estimate))%>%
      mutate(b=-log10(padj)*sign(estimate))%>%
      dplyr::select(factor, clin_var, b)%>%
      pivot_wider(names_from = clin_var,values_from = b)%>%
      filter(!is.na(factor))%>%
      column_to_rownames('factor')%>%as.matrix()

    clin_df_tests_pmat <- clin_df_tests_pmat[paste0('Factor',
                                                  sort(as.integer(gsub('Factor','', rownames(clin_df_tests_pmat))))),]
    clin_df_tests_emat <- clin_df_tests_emat[paste0('Factor',
                                                    sort(as.integer(gsub('Factor','', rownames(clin_df_tests_emat))))),]


    #Those with at least one significant hit
    onesig <- apply(clin_df_tests_pmat, 2, function(x) sum(x<0.05)) %>% .[.>0] %>% na.omit()

    pdf(paste0(outdir,'Factors_clin_test_reduced',append,'.pdf'),width = 12,height = 10)
    corrplot::corrplot(clin_df_tests_emat[,names(onesig)],
                       col = colorRampPalette(colors = c('Blue', 'white','darkred'))(100),
                       tl.col = 'black', is.corr = F, addgrid.col = NA, outline = F,mar=c(0,0,2,0),
                       p.mat = clin_df_tests_pmat[,names(onesig)],sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
                       insig = 'label_sig', pch.col = 'black')
    dev.off()

#  }, mc.cores=2)
#}, mc.cores=4)
