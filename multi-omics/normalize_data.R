suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
  library(reticulate)
  library(clusterProfiler)
  INT <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
  
})
covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"

#1. Process methylation to a file
meth1 <- readRDS('2000HIV/AnalysisDV/Validation_MethylationDataBased/1qc6/2000HIV.Mvalue.rds')
meth2 <- readRDS('2000HIV/AnalysisDV/Discovery_MethylationDataBased/1qc6/2000HIV.Mvalue.rds')

meth_int <- intersect(rownames(meth1),rownames(meth2))

meth1 <- meth1[meth_int,]
meth2 <- meth2[meth_int,]

meth <- cbind(meth1, meth2)

#Remove three duplicated IDs, Xun tells me that they remove the first three ones
meth <- meth[,!duplicated(colnames(meth), fromLast=T)]

#Only independent variables
#methcor <- findCorrelation(cor(t(meth)), cutoff = .7)

#Take the top 20k with highest variation. I'll use variance instead of CV cause I assume all values are scaled
meth_var <- apply(meth, 1, var)
meth_top_var <- names(sort(meth_var, decreasing = T)[1:20000])
meth <- meth[meth_top_var,]
#meth <- melt(meth, variable.name='feature') %>% mutate(view = 'meth')
meth_file <- paste0(outdir, 'meth.csv')
#Write and then read
fwrite(as.data.frame(meth), meth_file, row.names = T)


#2. Read all files
cytfile <- 'meta_cQTL/data/2000HIV/Phenotype/cytokines_modified.tsv'
metabfile <- 'meta_cQTL/data/2000HIV-EU-discovery/metabolites/phenotype.tsv'
protfile <-'meta_cQTL/data/2000HIV-EU-discovery/proteins/phenotype.tsv'
gexfile <- 'meta_cQTL/data/2000HIV-EU-discovery/expression/phenotype.tsv'

plot_dist <- function(file, norm = NULL, group = F){
  tmpdf <- fread(file) %>% column_to_rownames(colnames(.)[1]) %>% t()%>% as.data.frame()
  if(!is.null(norm)) tmpdf <- apply(tmpdf, 2, norm) %>% as.data.frame()
  if(group){
    ggplot(reshape2::melt(tmpdf),aes(x=value, group=variable))+geom_density()+
      theme_bw()+theme(panel.grid = element_blank(), strip.background = element_blank())
  }else{
    ggplot(reshape2::melt(tmpdf),aes(x=value))+geom_density()+
      theme_bw()+theme(panel.grid = element_blank(), strip.background = element_blank())
  }
}

plot_dist(cytfile, log1p)
plot_dist(metabfile, log1p)
plot_dist(protfile)
plot_dist(gexfile, log1p)
plot_dist(meth_file)

