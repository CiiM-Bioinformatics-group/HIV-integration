suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'

table <- fread(covs)  %>% filter(ETHNICITY.White == 1)
table_grouped <- table %>% 
  group_by(Cohort) %>% summarise_at(c('AGE','BMI_BASELINE','CD4_LATEST','HIV_DURATION'), function(x){
    x = na.omit(x)
    paste0(round(median(x),digits = 1),' [', round(min(x),1),', ', round(max(x),1),']')
  }) %>% column_to_rownames('Cohort') %>% t()

rownames(table_grouped) <- str_to_sentence(rownames(table_grouped)) %>% gsub('_', ' ', .) %>% gsub('Bmi','BMI', .)%>%
  gsub('Cd4','CD4',.) %>% gsub('Hiv', 'HIV',.)
  
table(table$SEX_BIRTH, table$Cohort)
