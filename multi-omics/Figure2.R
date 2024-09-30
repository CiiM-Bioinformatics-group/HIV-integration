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
  library(patchwork)
  library(corrplot)
  library(rstatix)
  library(lavaan)
})

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'figure2'))
append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Colors: From jbgb13/peRReo
source('2000HIV/cQTL/mofa/code/peRReo.R')
viewcol <- latin_palette('buenavista', n=5)
names(viewcol) <- views_names(model)

#The genetic component
mofaqtl <- fread('meta_cQTL/out/2000HIV-EU-discovery/mofa/mapping/main_1e-5.tsv')%>%
  group_by(gene) %>%slice_min(`p-value`,n = 1)

pthresh <- -log10(5e-8/length(unique(mofaqtl$gene)))

barplot_mofaqtl1 <- mofaqtl %>%
  mutate(gene = factor(gene, levels = paste0('Factor',length(unique(mofaqtl$gene)):1)))%>%
  ggplot(aes(y=gene, x=-log10(`p-value`),fill=-log10(`p-value`)>pthresh))+
  geom_col()+
  geom_text(aes(y=gene,x=1,label=gene,color=-log10(`p-value`)>pthresh))+
  theme_bw()+
  scale_fill_manual(values = c('ivory2', '#6A4D74'))+
  scale_color_manual(values = c('grey','white'))+
  #facet_grid(.~-log10(`p-value`)>10, scales = 'free')+
  geom_vline(aes(xintercept=pthresh), linetype = 'dashed')+
  xlab('-log10(P)')+
  scale_x_continuous(position = 'top')+
  coord_cartesian(xlim = c(0,10))+
  scale_y_discrete(position = 'left')+
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(),axis.text.y.left = element_blank(),
        axis.title.y = element_blank(), panel.border = element_blank(),legend.position = 0,strip.placement = 'outside',
        panel.spacing = unit(0, "lines"), strip.background = element_blank(),
        axis.line.x = element_line(),plot.margin = unit(c(0,0,0,0), 'cm'))
maxp <- max(-log10(mofaqtl$`p-value`))
barplot_mofaqtl <- barplot_mofaqtl1+
  plot_spacer()+
  barplot_mofaqtl1+
  coord_cartesian(xlim = c(maxp-2.5,maxp))+
  scale_x_continuous(position = 'top', breaks = c(round(maxp-2.5, digits = 1),round(maxp, digits = 1)))+
  xlab('')+
  plot_layout(widths = c(160,1,40))

ggsave(paste0(outdir,'figure2/genetics_barplot',append,'.pdf'),barplot_mofaqtl, width = 4, height = 4)


topsnp <- 'chr19:53824059:C:A;rs34436714'
topchr <- str_split(topsnp, pattern = ':', simplify=T)[,1]

xqtls <- fread(cmd=paste0('grep -w "',topsnp,
                      '" meta_cQTL/out/2000HIV-EU-discovery/*/mapping/main_genomewide.tsv'))%>%
  mutate(V1 = gsub('meta_cQTL/out/2000HIV-EU-discovery/','', V1)%>%
           gsub(paste0('/mapping/main_genomewide.tsv:',topsnp), '', .))%>%
  filter(V1 %in% c('proteins','metabolites','cell_count','expression'))%>%
  mutate(V1=gsub('_many','',V1))%>%
  mutate('sig'= 'P<5e-8')

#Replication of the effects
colnames(xqtls) <- c('view','gene','beta','tstat','pval','sig')
matches <- pull(xqtls, gene, name = view)
matches_list <- split(matches, names(matches))

xqtl_rep <- lapply(names(matches_list), function(v){
  print(v)
  g <- matches_list[[v]]

  emptydf <- data.frame('view' = v, 'gene' = g, 'beta' = NA, 'tstat' = NA, 'pval' = NA)
  file <- paste0('meta_cQTL/out/2000HIV-EU-validation/',v,'/mapping/main/',topchr,'.tsv')

  if(!file.exists(file)) {
    print('File does not exist')
    return(emptydf)
    }
  command <- paste0('grep ',paste(paste0('-e "',topsnp, '.', g,'"'), collapse = ' '),
                    ' meta_cQTL/out/2000HIV-EU-validation/',v,'/mapping/main/',topchr,'.tsv')
  df <- fread(cmd=command,
              col.names = c('V1',colnames(xqtls)[2:5]))
  if(is.null(df)) df <- emptydf else df <- dplyr::select(df,-V1)%>%mutate(view=v)

  return(df)
}) %>% bind_rows()

xqtl_with_rep <- left_join(xqtls, xqtl_rep, by=c('view','gene'), suffix = c('','.validation'))%>%
  mutate(fdr.validation = p.adjust(pval.validation, method = "fdr")) %>%
  mutate(replicated = (sign(beta) == sign(beta.validation) & fdr.validation < 0.05))
fwrite(xqtl_with_rep, paste0(outdir,'figure2/nlrp12_qtls_replicated',append,'.csv'))
xqtls_dotplot <- xqtl_with_rep %>% filter(replicated) %>%
  group_by(sig,view)%>%summarise(n=n())%>%
  mutate(view = ifelse(view == 'expression','gex',
                     ifelse(view == 'proteins', 'prot', ifelse(view == 'metabolites', 'metab', view))))%>%
  ggplot(aes(x=sig,y=n,color=view,fill=view))+
  geom_col(width=.03, position = position_dodge(width = .5))+
  geom_point(aes(size=n), position = position_dodge(width = .5))+
  geom_text(aes(label=n,y=n+10), position = position_dodge(width = .5))+
  scale_color_manual(values=c(viewcol, 'cell_count' = 'orange'))+
  scale_fill_manual(values=c(viewcol, 'cell_count' = 'orange'))+
  #facet_grid(rows='view')+
    theme_bw()+
  xlab('')+ylab(paste('xQTL for',topsnp))+
  theme(panel.grid = element_blank(),strip.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(),
        axis.line.x = element_line(),legend.position="top")+
  guides(size = "none")
ggsave(paste0(outdir,'figure2/nlrp12_qtls_dotplot_replicated',append,'.pdf'),xqtls_dotplot, width = 2, height = 5)

#The best inflammasome score
inflam_gmt <- read.gmt('inflammasome_gmt/all.gmt')
inflam_vec <- pull(inflam_gmt, gene, name=term)
inflam_list <- split(inflam_vec, names(inflam_vec))

#Calculate for gex and proteins
features_names(model)$prot <- features_names(model)$prot %>% gsub('_Inflam','-Inflam', .) %>% gsub('_Cardio','-Cardio',.)%>%
  gsub('_Onco','-Onco',.)%>%gsub('_Neuro','-Neuro',.)%>%gsub('_II$','-II',.) %>%
  str_split(pattern = '-', simplify = T) %>% .[,1]

gex_data <- get_data(model, views = 'gex')[[1]][[1]]
prot_data <- get_data(model, views = 'prot')[[1]][[1]]

gex_inflam <- gex_data[rownames(gex_data) %in% inflam_gmt$gene,]
prot_inflam <- prot_data[rownames(prot_data) %in% inflam_gmt$gene,]

gex_prot_samples <- intersect(colnames(gex_inflam), colnames(prot_inflam))

gex_int <- t(gex_inflam[,gex_prot_samples])
prot_int <- t(prot_inflam[,gex_prot_samples])

colnames(gex_int) <- paste0(colnames(gex_int), '_gex')
colnames(prot_int) <- paste0(colnames(prot_int), '_prot')

gex_prot_cor <- psych::cor2(gex_int,prot_int)

gex_clust <- hclust(dist(gex_prot_cor)) %>% cutree(3) #First one is  inflammasome, second one NFKB
gex_clust <- ifelse(gex_clust == 2, 'NFKB', ifelse(gex_clust == 1, 'Inflam', 'Other'))
prot_clust <- hclust(dist(t(gex_prot_cor))) %>% cutree(2) #First one is  inflammasome
prot_clust <- ifelse(prot_clust == 1, 'Inflam', 'Other')


#How do they correlate in general
pdf(paste0(outdir,'figure2/inflammasome_gex_to_prot_heatmap',append,'.pdf'), width = 10, height = 10)
print(Heatmap(gex_prot_cor, left_annotation = rowAnnotation('cl' = gex_clust),
        top_annotation = HeatmapAnnotation('cl' = prot_clust),
        show_row_dend = F, show_column_dend = F))
dev.off()

#A PCa on only these genes/proteins
#pc_inflam <- prcomp(cbind(gex_int,prot_int))
#Heatmap(pc_inflam$rotation[,1],left_annotation = rowAnnotation('cl' = c(gex_clust,prot_clust)),
#        cluster_rows = T)
#pc_score <- pc_inflam$x[,1] %>% as.data.frame() %>% rownames_to_column('sample')%>%
#  dplyr::rename('pc_score' = colnames(.)[2])


#inflam_list <- append(inflam_list, list('prot_inflam' = names(prot_clust[prot_clust == 'Inflam'])%>%gsub('_prot','',.),
#                    'gex_nfkb' = names(gex_clust[gex_clust == 'NFKB'])%>%gsub('_gex','',.),
#                    'gex_inflam' = names(gex_clust[gex_clust == 'Inflam'])%>%gsub('_gex','',.),
#                    'gex_all' = names(gex_clust)%>%gsub('_gex','',.)))

inflam_scores <-  lapply(names(inflam_list), function(gl){
  print(gl)
    if(grepl('prot',gl)){
      score <- apply(prot_inflam[inflam_list[[gl]],], 2, mean)
      return(data.frame('sample' = names(score),'score' = score, 'n' = length(inflam_list[[gl]]), 'path' = gl))
    }else if(grepl('gex',gl)){
      score <- apply(gex_inflam[inflam_list[[gl]],], 2, mean)
      return(data.frame('sample' = names(score),'score' = score, 'n' = length(inflam_list[[gl]]), 'path' = gl))
    }else{
      genes_avail <- intersect(inflam_list[[gl]], rownames(gex_inflam))
      gex_score <- apply(gex_inflam[genes_avail,], 2, mean)
      prots_avail <- intersect(inflam_list[[gl]], rownames(prot_inflam))
      if(length(prots_avail) < 2) return(data.frame('sample' = NA,'score' = NA, 'n' = NA, 'path' = paste0('gex_',gl)))
      prot_score <- apply(prot_inflam[prots_avail,], 2, mean)
      df <- rbind(data.frame('score' = gex_score, 'n' = length(genes_avail),'path' = paste0('gex_',gl)) %>%
                    rownames_to_column('sample'),
                  data.frame('score' = prot_score,'n' = length(prots_avail),'path' = paste0('prot_',gl)) %>%
                    rownames_to_column('sample'))
    }
  })%>%bind_rows()%>%
  separate(path, into = c('view', 'path'), extra = 'merge', sep = '_')#%>%  left_join(pc_score)

fwrite(inflam_scores,paste0(outdir,'figure2/inflammasome_scores',append,'.csv'))

#Correlate both
ggplot(inflam_scores %>% dplyr::select(-n)%>%
         pivot_wider(names_from = view, values_from = score)%>%
         filter(!is.na(prot)), aes(x=prot,y=gex))+
  geom_point(alpha=.2)+
  geom_smooth(method = 'lm')+
  scale_color_viridis_d()+
  stat_cor()+
  facet_wrap(~path)+
  theme_bw()+
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())
ggsave(paste0(outdir,'figure2/inflammasome_gex_to_prot',append,'.pdf'), width = 4, height = 4)

#All with PCs
#ggplot(inflam_scores,
#       aes(x=score,y=pc_score,color=path))+
#  geom_point(alpha=.2)+
#  geom_smooth(method = 'lm')+
#  facet_wrap(~view, scales = 'free')+
#  stat_cor()+
#  theme_bw()+
#  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line())

#Correlate both to Factor 16
inflam_scores_factor <- inflam_scores %>%
  left_join(get_factors(model,as.data.frame=T))

topfactor <- filter(mofaqtl, `p-value` == min(mofaqtl$`p-value`)) %>% pull(gene)

ggplot(inflam_scores_factor %>% filter(factor == topfactor), aes(x=value,y=score))+
  geom_point(alpha=.2)+
  geom_smooth(method = 'lm')+
  stat_cor()+xlab(topfactor)+
  facet_grid(view~path, scales = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())
ggsave(paste0(outdir,'figure2/inflammasome_score_to_',topfactor,append,'.pdf'), width = 8, height = 4)

ggplot(inflam_scores_factor %>% filter(factor == topfactor & path == 'REACTOME_THE_NLRP3_INFLAMMASOME'),
       aes(x=value,y=score))+
  geom_point(alpha=.2)+
  geom_smooth(method = 'lm')+
  stat_cor()+xlab(topfactor)+
  facet_grid(view~path, scales = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())
ggsave(paste0(outdir,'figure2/inflammasome_score_reactome_to_',topfactor,append,'.pdf'), width = 3.5, height = 6)

ggplot(inflam_scores_factor %>% filter(path == 'REACTOME_THE_NLRP3_INFLAMMASOME'),
       aes(x=value,y=score,color=view,path))+
  geom_point(alpha=.05)+
  geom_smooth(method = 'lm')+
  scale_color_manual(values = viewcol)+
  stat_cor()+
  facet_wrap(~factor)+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(paste0(outdir,'figure2/inflammasome_score_reactome_to_allfactors',append,'.pdf'), width = 12, height = 12)

#Check if the qtl snp has an effect on both
topsnp2 <- gsub('chr','',topsnp) %>% str_split(pattern = ';', simplify=T)%>%.[,1]
snp_dos <- fread(cmd = paste0('zgrep -w -e CHROM -e "', topsnp2,
                               '" meta_cQTL/data/2000HIV-EU-discovery/Genotype/vcf/',topchr,'.vcf.gz'))%>%
  as.data.frame()
snp_dos <- snp_dos[,10:ncol(snp_dos)]
colnames(snp_dos) <- str_split(colnames(snp_dos), pattern = '_', simplify = T)[,2]
snp_dos <- mutate_all(snp_dos, function(x) str_split(x, pattern = ':', simplify = T)[,2] %>% as.numeric())
snp_dos <- t(snp_dos) %>% as.data.frame() %>% rownames_to_column('sample')
snp_dos <- dplyr::rename(snp_dos, 'SNP' = V1)

topalleles <- str_split(topsnp2,pattern = ':', simplify = T)[,c(3,4)]

inflam_scores_factor_snp <- left_join(inflam_scores_factor, snp_dos, by = 'sample')%>%
  mutate(snp = topsnp)%>%separate(snp,into = c('chr','pos','ref','alt'))%>%
  mutate(GT = ifelse(SNP > 1.5, paste0(alt,alt), ifelse(SNP > 0.5, paste0(ref,alt), paste0(ref,ref))))%>%
  filter(!is.na(SNP))%>%
  mutate(GT = factor(GT, levels = c(paste0(topalleles[1],topalleles[1]),paste0(topalleles[1],topalleles[2]),
                                    paste0(topalleles[2],topalleles[2]))))

ggplot(inflam_scores_factor_snp %>%
         dplyr::select(SNP,score,view,path) %>% distinct(),
       aes(x=SNP, y = score))+
  geom_point()+
  geom_smooth(method = 'lm')+
  stat_cor()+
  facet_grid(view~path, scales = 'free_y')+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(paste0(outdir,'figure2/inflammasome_score_to_snp_dosage',append,'.pdf'), width = 8, height = 4)


ggboxplot(inflam_scores_factor_snp %>%
            dplyr::select(GT,score,view,path) %>% distinct(),
          x='GT', y = 'score', add = 'jitter', add.params = list(alpha=.05))+
  stat_pwc(label = 'p.signif', hide.ns = T)+
  facet_grid(view~path, scales = 'free_y')
ggsave(paste0(outdir,'figure2/inflammasome_score_to_snp_gt',append,'.pdf'), width = 8, height = 4)

#ggboxplot(inflam_scores_factor_snp %>%
#            dplyr::select(GT,pc_score) %>% distinct(),
#          x='GT', y = 'pc_score', add = 'jitter', add.params = list(alpha=.05))+
#  stat_pwc(label = 'p.signif', hide.ns = T)

ggboxplot(inflam_scores_factor_snp %>%
  filter(path == 'REACTOME_THE_NLRP3_INFLAMMASOME')%>%
            dplyr::select(GT,score,view,path) %>% distinct(),
          x='GT', y = 'score', add = 'jitter', add.params = list(alpha=.05))+
  stat_pwc(label = 'p.signif', hide.ns = T)+
  facet_grid(view~path, scales = 'free_y')+
  xlab(topsnp)+theme(strip.background = element_blank())
ggsave(paste0(outdir,'figure2/inflammasome_score_reactome_to_snp_gt',append,'.pdf'), width = 3, height = 4)


ggplot(inflam_scores_factor_snp %>%
         filter(path == 'REACTOME_THE_NLRP3_INFLAMMASOME')%>%
         dplyr::select(SNP,value,factor) %>% distinct(),
       aes(x=SNP, y = value))+
  geom_point(alpha=.05)+
  geom_smooth(method = 'lm')+
  stat_cor()+
  facet_wrap(~factor, scales = 'free_y')+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(paste0(outdir,'figure2/allfactors_to_snp_dosage',append,'.pdf'), width = 12, height = 12)

#Correlate both to the metabolites we have found
metab_ids  <- fread('2000HIV/Metabolites/ionMz2ID.tsv')%>%
  filter(formula %in% (xqtls %>% filter(view == 'metabolites' & pval<5e-8) %>% pull(gene)))

inflam_metabs <- get_data(model, views = 'metab', as.data.frame=T)%>%
  filter(feature %in% metab_ids$name)

inflam_scores_metab <- inner_join(inflam_scores, inflam_metabs %>% dplyr::select(-view), by = 'sample')

ggplot(inflam_scores_metab%>%
         dplyr::select(score,value,view,path,feature) %>% distinct(),
       aes(x=score,y=value, color=feature))+
  geom_point(alpha = .05)+
  geom_smooth(method = 'lm')+
  stat_cor()+
  scale_color_manual(values = c('darkred','darkgreen'))+
  facet_grid(view~path, scales = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())
ggsave(paste0(outdir,'figure2/inflammasome_score_to_metabolites_nlrp12',append,'.pdf'), width = 8, height = 4)

ggplot(inflam_scores_metab%>%
filter(path == 'REACTOME_THE_NLRP3_INFLAMMASOME')%>%
         dplyr::select(score,value,view,path,feature) %>% distinct(),
       aes(y=score,x=value, color=feature))+
  geom_point(alpha = .05)+
  geom_smooth(method = 'lm')+
  stat_cor()+
  scale_color_manual(values = c('black','darkgrey','lightgrey'))+
  facet_grid(view~path, scales = 'free')+
  theme_bw()+
  theme(panel.grid = element_blank(), strip.background = element_blank())
ggsave(paste0(outdir,'figure2/inflammasome_score_reactome_to_metabolites_nlrp12',append,'.pdf'), width = 4, height = 4)

#Correlate the factor to cell counts
cellcounts <- fread('meta_cQTL/data/2000HIV-EU-discovery/cell_count/phenotype.tsv')%>%
  column_to_rownames('V1') %>% .[filter(xqtls, view == 'cell_count') %>% pull(gene),]%>%
  t()%>%as.data.frame()%>%rownames_to_column('sample')%>%melt(value.name = 'cc')

#Are percentages relative to all mono? Should be cause all are almost 100
fread('meta_cQTL/data/2000HIV-EU-discovery/cell_count/phenotype.tsv')%>%
  column_to_rownames('V1') %>% .[filter(xqtls, view == 'cell_count') %>% pull(gene),]%>%t()%>%rowSums()

cellcounts_factors <- inner_join(cellcounts,
                                get_factors(model,factors = as.numeric(gsub('Factor','',topfactor)),
                                            as.data.frame=T), by='sample')

ggplot(cellcounts_factors, aes(x=value,y=cc))+
  geom_point(alpha=.1)+geom_smooth(method = 'lm')+
  stat_cor()+
  facet_wrap(~variable)+theme_bw()+theme(panel.grid = element_blank(), strip.background = element_blank())+
  xlab(topfactor)+ylab('Cell counts')
ggsave(paste0(outdir,'figure2/cell_counts_to_',topfactor,append,'.pdf'), width = 6, height = 3)

#Is the inflammasome score different in people with PLAQUE?
clin <- fread(paste0(outdir,'clin_cor/clin_df_corrected_scaled.csv')) %>%
  dplyr::rename('sample' = colnames(.)[1])

inflam_scores_plaque <- inner_join(inflam_scores, clin %>% dplyr::select(sample,PLAQUE_YESNO),by='sample')%>%
  filter(!is.na(PLAQUE_YESNO))

ggboxplot(inflam_scores_plaque, x = 'PLAQUE_YESNO', y = 'score')+
  facet_grid(view~path, scales = 'free_y')+stat_pwc(hide.ns = T, label = 'p.signif')+
  theme(strip.background = element_blank())
ggsave(paste0(outdir,'figure2/inflammasomescore_to_plaque.pdf'), width = 6, height = 3)

#Mediation for snp to gene expression to proteins
# mediation_df <- inflam_scores_factor_snp %>%
#   filter(path == 'REACTOME_THE_NLRP3_INFLAMMASOME')%>%
#   dplyr::select(sample,SNP,score,view)%>%
#   distinct()%>%
#   pivot_wider(names_from = 'view', values_from = 'score')%>%
#   column_to_rownames('sample')
#
# cor(mediation_df)
#
# library(mediation)
# results_prot <- mediate(
#   lm(prot ~ SNP, mediation_df),
#   lm(gex ~ SNP+prot, mediation_df),
#   treat = 'SNP',
#   mediator = 'prot',
#   boot = TRUE,
#   sims = 500
# )
# results_gex <- mediate(
#   lm(gex ~ SNP, mediation_df),
#   lm(prot ~ SNP+gex, mediation_df),
#   treat = 'SNP',
#   mediator = 'gex',
#   boot = TRUE,
#   sims = 500
# )
# summary(results_prot)
# summary(results_gex)
#
# #What is the R2 of a model in which we predict il1b with age, sex, genetics (X variants from cQTL)
# #and mofa factors (including Factor6). Also in validation
