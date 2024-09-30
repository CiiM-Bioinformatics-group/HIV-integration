#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(coloc)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)


dataroot <- "meta_cQTL/data"
root <- "meta_cQTL/out"
cohort <- "2000HIV-EU-discovery"
cov <- "main"
pheno <- "protein"

outdir <- paste0(root, "/", cohort, "/", pheno, "/extra")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pval <- "studywide"

# read in our qtl
loci <- TRUE
if (loci) {
    qtl <- fread(paste0(root, "/", cohort, "/", pheno, "/mapping/", cov, "_", pval, "_", "loci.tsv")) #%>%
    #mutate(LOG10P = -log10(pval))
} else {
    qtl <- fread(paste0(root, "/", cohort, "/", pheno, "/mapping/", cov, "_1e-5.tsv")) %>%
    separate("SNP", into = c("CHR", "BP", "REF", "ALT", "rsid"), remove = TRUE) #%>%
    #mutate(LOG10P = -log10(pval))
}
nrow(qtl) #

qtl_other <- fread("resources/GTEx/v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")






#OLD
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
})

outdir <- '2000HIV/cQTL/mofa/out/pqtl_ukb/'
dir.create(outdir)

pqtl <- fread('meta_cQTL/out/2000HIV-EU-discovery/proteins/extra/main_genomewide_loci_ukb.tsv')%>%
  separate(ID, into = c('chr_hg19','pos_hg19','REF_ukb','ALT_ukb'), sep = ':', remove = F)%>%
  mutate('beta_ukb_aligned' = ifelse(REF == REF_ukb & ALT == ALT_ukb, beta_ukb,
                                     ifelse(REF == ALT_ukb & ALT == REF_ukb, -beta_ukb,NA)),
         'se' = beta/`t-stat`)%>%
  mutate(prot = gene %>% gsub('_Inflam','-Inflam', .) %>% gsub('_Cardio','-Cardio',.)%>%
           gsub('_Onco','-Onco',.)%>%gsub('_Neuro','-Neuro',.)%>%gsub('_II$','-II',.))%>%
  separate(prot, into = c('prot','type'), sep = '-')%>%
  mutate('t_diff_effect' = (beta-beta_ukb)/(sqrt(se^2+se_ukb^2)),
         'p_diff_effect' = 2 * (pt(abs(t_diff_effect), 1051, lower.tail = F)),
         'padj_diff_effect' = p.adjust(p_diff_effect, method = 'fdr'))


ggplot(pqtl, aes(x=beta, y=beta_ukb_aligned))+
  geom_point()+
  geom_smooth(method = 'lm')+
  geom_vline(aes(xintercept=0), linetype='dashed')+
  geom_hline(aes(yintercept=0), linetype='dashed')+
  stat_cor()+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(paste0(outdir,'pqtl_beta_cor.pdf'),width = 4,height = 4)

pqtl_discordant <- pqtl %>% filter(sign(beta) == 1 & sign(beta_ukb_aligned) == -1) %>%
  filter(pval_ukb < 1e-8)
fwrite(pqtl_discordant, paste0(outdir,'pqtl_discordant.csv'))

pqtl %>%
  ggplot(aes(x=t_diff_effect, fill = type))+
  geom_density()+
  facet_wrap(~type, ncol = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(paste0(outdir,'pqtl_t_density.pdf'),width = 6,height = 4)

pqtl %>%
  ggplot(aes(x=t_diff_effect, fill = type))+
  geom_density()+
  facet_grid(sign(`t-stat`)~type)+
  theme_bw()+
  theme(panel.grid = element_blank())
pqtl %>%
  ggplot(aes(y=abs(t_diff_effect), x=abs(beta)))+
  geom_point()+
  geom_smooth(method = 'lm')+
  stat_cor()+
  facet_grid(~(locus==1974))+
  theme_bw()+
  theme(panel.grid = element_blank())

pqtl %>%
  filter(locus != 1974)%>%
  ggplot(aes(x=t_diff_effect, fill = type))+
  geom_density()+
  facet_wrap(~type, ncol = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(paste0(outdir,'pqtl_t_density_no_nlrp12.pdf'),width = 6,height = 4)

qqman::qq(pqtl$padj_diff_effect)
qqman::qq(filter(pqtl, locus != 1974) %>% pull(padj_diff_effect))

pqtl %>%
  mutate(t_ukb = beta_ukb/se_ukb, MAF = ifelse(A1FREQ > 0.5, 1-A1FREQ,A1FREQ))%>%
  rowwise()%>%
  mutate(min_t = min(abs(t_ukb),abs(`t-stat`)))%>%
  ggplot(aes(x=MAF, y=abs(t_ukb)))+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(paste0(outdir,'pqtl_maf.pdf'),width = 4,height = 4)


ci <- qt(0.975,1000)

pqtl %>%
  filter(locus == 1974)%>%
  dplyr::select(beta,se,beta_ukb,se_ukb,prot)%>%
  rename(beta_hiv=beta,se_hiv=se)%>%
  pivot_longer(c(beta_hiv,beta_ukb, se_hiv, se_ukb), names_to = c('.value','cohort'),
               names_pattern = '([a-z]+)_([a-z]+)')%>%
  mutate(ci_low = beta-ci*se, ci_high=beta+ci*se)%>%
  ggplot()+
  geom_point(aes(x=beta,color=cohort,y=prot))
ggsave(paste0(outdir,'pqtl_beta_nlrp12.pdf'),width = 4,height = 8)


pqtl_perlocus <- pqtl %>%
  group_by(locus) %>%
  summarise(median_t = median(t_diff_effect), n=n(), median_beta = median(beta))

ggplot(pqtl_perlocus,aes(x=locus,y=median_t))+
  geom_col(fill='darkred')+
  theme_bw()+theme(panel.grid = element_blank())
ggsave(paste0(outdir,'pqtl_locus_median_t.pdf'),width = 5,height = 5)

ggplot(pqtl_perlocus,aes(x=median_t))+
  geom_density()+
  theme_bw()+theme(panel.grid = element_blank())
ggsave(paste0(outdir,'pqtl_locus_median_t_density.pdf'),width = 4,height = 4)
ggplot(pqtl_perlocus,aes(x=abs(median_beta),y=abs(median_t)))+
  geom_point()+
  geom_smooth(method = 'lm')+
  stat_cor()+
  theme_bw()+theme(panel.grid = element_blank())
