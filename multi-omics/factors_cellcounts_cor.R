suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(ComplexHeatmap)
  library(ggpubr)
  library(parallel)
  library(rstatix)
  library(xlsx)
})

covs = '/vol/projects/CIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "/vol/projects/CIIM/2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'cc_cor'))
#mclapply(c('corrected','uncorrected'), function(i){
#  mclapply(c('scaled','unscaled'), function(j){
    #append = paste0('_',i,'_',j)
    append = '_corrected_scaled'
    #args = commandArgs(trailingOnly=TRUE)
    #if(is.na(args[1])) corrected = '' else corrected = args[1]
    #modelfile = paste0(outdir,"model_corrected_unscaled.hdf5")

    #Load
    model <- readRDS(paste0(outdir,'model', append, '.rds'))

    #Factors
    factors <- get_factors(model, factors = "all") %>% .$single_group

cc <- readxl::read_excel('/vol/projects/CIIM/2000HIV/2000HIV_Flowcytometry/UpdatedFlowcytometry/2000HIV_FLOW_ABS_panel123merged_QCed_untransformed(raw)data_1423samples_356vars_Nov062023.xlsx',
                          1)

factors_cc <- inner_join(as.data.frame(factors) %>% rownames_to_column('samples'),
                             rename(cc, 'samples' = SampleID))%>%
      column_to_rownames('samples')



factors_cc_tests <- factors_cc %>%
      pivot_longer(colnames(.)[grepl('Factor',colnames(.))], names_to = 'factor', values_to = 'factor_value')%>%
      pivot_longer(colnames(.)[!grepl('factor', colnames(.))], names_to = 'cc', values_to = 'cc_value')%>%
      group_by(factor,cc)%>%
      cor_test(vars1 = 'factor_value', vars2 = 'cc_value', method = 'spearman')

factors_cc_tests$padj <- p.adjust(factors_cc_tests$p, method = 'bonferroni')

factors_cc_tests <- factors_cc_tests %>% arrange(padj)

fwrite(factors_cc_tests, paste0(outdir,'cc_cor/factors_cc_test',append,'.csv'), row.names = T)

    #Single correlation
factors_cc_cor <- psych::corr.test(
      factors_cc[,grepl('^Factor',colnames(factors_cc))],
      factors_cc[,!grepl('^Factor',colnames(factors_cc))],
      use = 'pairwise.complete')

#fwrite(as.data.frame(factors_cov_cor[c('r','p','p.adj')]), paste0(outdir,'factors_cov_cor',append,'.csv'), row.names = T)

#Only the significant ones
cors <- t(factors_cc_cor$r)
sigcors <- factors_cc_tests %>% dplyr::select(factor,cc,padj)%>%
      pivot_wider(names_from = factor, values_from=padj)%>%column_to_rownames('cc')

sigcors <- sigcors[rownames(cors),colnames(cors)]
cors[sigcors > 0.001] <- 0
cors <- cors[!rowSums(cors) == 0, !colSums(cors) == 0]

pdf(paste0(outdir,'cc_cor/Factors_cc_cor',append,'.pdf'),width = 10,height = 40)
Heatmap(cors, show_row_dend = F, cluster_columns = F)
dev.off()

#Only the CCR5 ones
pdf(paste0(outdir,'cc_cor/Factors_cc_cor_ccr5',append,'.pdf'),width = 10,height = 10)
Heatmap(cors[grepl('CCR5', rownames(cors)),], show_row_dend = F, cluster_columns = F)
dev.off()
#Only the CXCR4 ones
pdf(paste0(outdir,'cc_cor/Factors_cc_cor_cxcr4',append,'.pdf'),width = 10,height = 10)
Heatmap(cors[grepl('CXCR4', rownames(cors)),], show_row_dend = F, cluster_columns = F)
dev.off()

#  }, mc.cores=2)
#}, mc.cores=4)


#A CCR5+ to CCR5- score
cc_ccr5 <- factors_cc %>% .[,grepl('CCR5|^Factor',colnames(.))]   %>%
  mutate(
    cd4ratio = `Panel2_2096_CD4+CXCR4+CCR5+`/`Panel2_2095_CD4+CXCR4+CCR5-`,
    cd8ratio = `Panel2_2100_CD8+CXCR4+CCR5+`/`Panel2_2099_CD8+CXCR4+CCR5-`
  )

cd4_test <- cc_ccr5 %>% dplyr::select(cd4ratio,colnames(.)[grepl('Factor',colnames(.))]) %>%
  rstatix::cor_test(cd4ratio)%>% mutate(padj = p.adjust(p, method = 'fdr'))

cd8_test <- cc_ccr5 %>% dplyr::select(cd8ratio,colnames(.)[grepl('Factor',colnames(.))]) %>%
  rstatix::cor_test(cd8ratio)%>% mutate(padj = p.adjust(p, method = 'fdr'))

ggplot(cc_ccr5, aes(x=cd4ratio,y=Factor16))+
  geom_point()+geom_smooth(method = 'lm')+stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())

ggplot(cc_ccr5, aes(x=`Panel2_2096_CD4+CXCR4+CCR5+`,y=Factor16))+
  geom_point()+geom_smooth(method = 'lm')+stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())

ggplot(cc_ccr5, aes(x=`Panel2_2095_CD4+CXCR4+CCR5-`,y=Factor16))+
  geom_point()+geom_smooth(method = 'lm')+stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())
