suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(ComplexHeatmap)
  library(ggpubr)
  library(rstatix)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
append = "_corrected_scaled"
alpha <- 0.01
append <- paste0('_validation_alpha',alpha,append)
outdir = "2000HIV/cQTL/mofa/out/"
#Read the clean table created by clin_df.R
clin_df_full <- fread(paste0(outdir,'validation/clin_df',append,'.csv'))%>%dplyr::rename('Record.Id' = colnames(.)[1])%>%
  left_join(fread(covs))%>%
      column_to_rownames(colnames(.)[1])

clin_df_clean <- clin_df_full %>% dplyr::select(colnames(.)[!colnames(.) %in% colnames(fread(covs, nrows = 1))])


    #Just correlation
    clin_df_cor <- psych::corr.test(
      clin_df_clean[,grepl('^Factor',colnames(clin_df_clean))],
      clin_df_clean[,!grepl('^Factor|^cluster',colnames(clin_df_clean))],
      use = 'pairwise.complete')

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
      obs <- length(unique(clin_df_clean[,var]))
      empty_df <- data.frame('p' = NA, 'factor' = NA, 'clin_var' = var, 'test' = NA)

      tmp.df <- clin_df_clean %>% dplyr::select(all_of(c(facts, var)))
      colnames(tmp.df)[ncol(tmp.df)] <- 'variable'
      
      tmp.df <- tmp.df %>% pivot_longer(cols = all_of(facts),names_to = 'factor', values_to = 'value')
      
      if(obs < 2){
        return(empty_df)
      }else if(obs== 2){
        out <- tmp.df %>%
          group_by(factor)%>%
          wilcox_test(value ~ variable, detailed = T, ref.group = '1')%>%
          dplyr::select(factor, p,estimate, group1, group2, n1, n2,conf.low,conf.high)%>%
          mutate('clin_var' = var, 'test' = 'wilcox')
      }else if(obs < 10){
        out <- tmp.df %>%
          group_by(factor)%>%
          kruskal_test(value ~ variable)%>%
          dplyr::select(factor, p)%>%mutate('clin_var' = var, 'test' = 'kruskal','estimate'=NA,'group1'=NA,
          'group2'=NA,'n1'=NA,'n2'=NA,'conf.low'=NA,'conf.high'=NA)
      }else{
        return(empty_df)
        }
    }) %>% bind_rows()

    #Correct
    clin_df_tests_mat <- na.omit(clin_df_tests) %>% mutate('padj' = p.adjust(p, method = 'holm'))%>%
      dplyr::select(factor, clin_var, padj)%>%
      pivot_wider(names_from = clin_var,values_from = padj)%>%
       column_to_rownames('factor')%>%as.matrix()

    clin_df_tests_mat <- clin_df_tests_mat[paste0('Factor',
                                                  sort(as.integer(gsub('Factor','', rownames(clin_df_tests_mat))))),]

    #Those with at least one significant hit
    onesig <- apply(clin_df_tests_mat, 2, function(x) sum(x<0.05)) %>% .[.>0] %>% na.omit()

if (length(onesig) > 0){
  pdf(paste0(outdir,'Factors_clin_test_reduced',append,'.pdf'),
      width = 10,height = 10)
  corrplot::corrplot(clin_df_tests_mat %>% apply(2,function(x) -log10(x))%>%.[,names(onesig)] %>% as.matrix(),
                     col = colorRampPalette(colors = c('white','darkred'))(100),
                     tl.col = 'black', is.corr = F, addgrid.col = NA, outline = F,mar=c(0,0,2,0),
                     p.mat = clin_df_tests_mat[,names(onesig)],sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,
                     insig = 'label_sig', pch.col = 'black')
dev.off()
}



    #Test specific examples:
    #Select only Nevirapine, INR_YES_NO, HIV stage

    colnames(clin_df_full) <- gsub('%', '_', colnames(clin_df_full))
    lapply(list(c('PLAQUE_YESNO', 'Factor6'),
                c('MH_TR_RESP_D_COPD', 'Factor11'),
                c('MH_TR_CV_D_Myocardial_infarction', 'Factor8'),
                c('RP_YESNO', 'Factor20')), function(z){
      clin_df_full %>%
      dplyr::select(c(z[1], z[2])) %>%
      na.omit() %>%
      ggboxplot(x=z[1], y= z[2], add = "jitter") +
       stat_pwc()
      ggsave(paste0(outdir,'validation/validation_boxplot_', paste(z[1:2], collapse='_'), append, '.pdf'),width=2,height=2)
      if(length(z)==3){
        clin_df_full %>% ggboxplot(x=z[1], y= z[2])+ stat_pwc()+facet_grid(c(z[3]))
        ggsave(paste0(outdir,'validation/validation_boxplot_', paste(z, collapse='_'), append, '.pdf'),width=4,height=4)
      }
    })
