suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'figure2'))
append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Read data
mofaqtl <- fread('meta_cQTL/out/2000HIV-EU-discovery/mofa/mapping/main_studywide.tsv')%>%
  group_by(gene) %>%slice_min(`p-value`,n = 1)

topsnp <- 'chr19:53824059:C:A;rs34436714'
topchr <- str_split(topsnp, pattern = ':', simplify=T)[,1]

xqtls <- fread(cmd=paste0('grep -w "',topsnp,
                          '" meta_cQTL/out/2000HIV-EU-discovery/*/mapping/main/',topchr,'.tsv'))%>%
  mutate(V1 = gsub('meta_cQTL/out/2000HIV-EU-discovery/','', V1)%>%
           gsub(paste0('/mapping/main/',topchr,'.tsv:',topsnp), '', .))

xqtls <- xqtls %>%
  filter(V1 %in% c('cytokines','proteins','metabolites','cell_count','expression_many'))%>%
  mutate(V1=gsub('_many','',V1))

colnames(xqtls) <- c('view','gene','beta','tstat','P')

xqtls %>% group_by(view) %>% summarise('Intercept' = summary(lm(tstat ~ 1))$coefficients[1,4],
                                       'Beta' = summary(lm(tstat ~ 1))$coefficients[1,1])

xqtls$padj <- p.adjust(xqtls$P)
