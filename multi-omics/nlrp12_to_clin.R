#Writing a table with the selected clinical variables and the factors calculated
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(parallel)
})

covs = 'meta_cQTL/data/2000HIV-EU-discovery/covariates.tsv'
outdir = "2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'clin_cor'))

clin_df_clean <- fread(paste0(outdir,'clin_cor/clin_df_corrected_scaled.csv'))%>%
  column_to_rownames(colnames(.)[1]) %>% .[,!grepl('^Factor',colnames(.))]%>%
  rownames_to_column('sample')
colnames(clin_df_clean) <- gsub('%','_', colnames(clin_df_clean))

nlrp12 <- fread(cmd = paste0('grep -e SNP -e "chr19:53824059:C:A;rs34436714"',
                ' meta_cQTL/out/2000HIV-EU-discovery/genotype/dosage/chr19.txt'))%>%
  column_to_rownames('SNP')%>%t()%>%as.data.frame()%>%dplyr::rename('SNP' = `chr19:53824059:C:A;rs34436714`)%>%
  rownames_to_column('sample')

covdf <- fread(covs, header = T) %>% column_to_rownames('V1') %>%t()%>%as.data.frame()%>%rownames_to_column('sample')

test <- inner_join(clin_df_clean, nlrp12)%>%inner_join(covdf, by = 'sample')


logistic <- lapply(colnames(clin_df_clean)[-1], function(p){
  formula <- as.formula(paste0(p,' ~ SNP + age + sex + bmi + lockdown + hiv_duration + season_cos + season_sin ',
                               '+ center_EMC + center_OLV'))
  print(formula)
  m <- glm(formula, data = test, family = 'binomial')
  df <- summary(m)$coefficients[2,]
  df <- as.data.frame(t(df))%>%mutate('clin' = p)
})%>%bind_rows()%>%mutate(padj = p.adjust(`Pr(>|z|)`, method = 'fdr'))

test <- mutate(test, GT = ifelse(SNP<0.5,'CC',ifelse(SNP<1.5,'CA','AA')) %>%
                 factor(levels = c('CC','CA','AA')))

table(test$GT,test$PLAQUE_YESNO)
chisq.test(table(test$GT,test$PLAQUE_YESNO))

test %>% dplyr::select(GT,PLAQUE_YESNO) %>% na.omit() %>%
  ggplot(aes(x=GT,fill=as.factor(PLAQUE_YESNO)))+geom_bar(position = 'fill')+
  theme_bw()+theme(panel.grid = element_blank(), panel.border = element_blank(),
                   axis.line.x = element_line(), axis.line.y = element_line())+
  scale_y_continuous(n.breaks = 3)+
  scale_fill_manual(values = c('grey','black'))+ylab('Proportion')+xlab('chr19:53824059:C:A;rs34436714')+
  labs(subtitle = paste0('Logistic regression p-value: ',
                         format(logistic$`Pr(>|z|)`[logistic$clin == 'PLAQUE_YESNO'],digits = 2)
                         ))
ggsave(paste0(outdir,'figure2/nlrp12_snp_to_plaque.pdf'), width = 5, height = 3)

logistic %>%
  filter(grepl('MH_TR_CV_', clin))%>%
  lm(`z value` ~ 1, data = .) %>%
  summary()
