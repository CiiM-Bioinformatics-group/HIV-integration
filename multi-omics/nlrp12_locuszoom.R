suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ggpubr)
  library(MOFA2)
  library(susieR)
  library(biomaRt)
  library(patchwork)
  library(ggrastr)
  library(ggrepel)
})

covs = '/vol/projects/CIIM/2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "/vol/projects/CIIM/2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'figure2'))
append = '_corrected_scaled'

#Load
model <- readRDS(paste0(outdir,'model', append, '.rds'))

#Read data
mofaqtl <- fread('/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/mofa/mapping/main_studywide.tsv')%>%
  group_by(gene) %>%slice_min(`p-value`,n = 1)

variant <- mofaqtl$SNP[which.min(mofaqtl$`p-value`)]
variant_split <- str_split(variant,pattern = '[:;]', simplify = T)
chr <- variant_split[,1]
pos <- as.numeric(variant_split[,2])
gene <- mofaqtl$gene[which.min(mofaqtl$`p-value`)]
dist = 5e5
min <- pos - dist
max <- pos + dist
path <- "/vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/mofa/mapping"
covar <- "main"

cmd <- paste0("awk -F ':'  '$2 > ", min, " && $2 < ", max, "' ", path, "/", covar, "/", chr, ".tsv", " | grep ", gene)
whole_locus <- fread(cmd = cmd, col.names = c("SNPid", "gene", "beta", "t", "P"))%>%
  separate(SNPid, into = c("CHR", "BP", "A1", "A2", "rsID"), sep = "[:;]", remove = F)


#Finemapping
#Extracting the dosage and correlation
dosage <- fread(cmd = paste0('grep "',whole_locus$SNPid[1],'" -A ',nrow(whole_locus)-1,
                       ' /vol/projects/CIIM/meta_cQTL/out/2000HIV-EU-discovery/genotype/dosage/',chr,'.txt'))
dosage <- as.data.frame(dosage)%>%column_to_rownames("V1")
R <- cor(t(dosage))
Z <- qnorm(whole_locus$P/2, lower.tail = F)*sign(whole_locus$beta)

#Plot
#p1 <- susie_plot(Z, y = "z", b=whole_locus$beta)

#Fine mapping
fitted_rss <- susie_rss(Z, R, L = 50,n = 1000)

#Transform into df
Susiedf <- data.frame("PIP" = fitted_rss$pip,  "CS" = NA)%>%
  rownames_to_column("SNP")%>%
  separate(SNP, into =c("CHR","BP","REF","ALT","rsID"),remove = F)


#Single locuszoom per gene
locuszoom <- function(point_color = "darkgrey", build, covar = 'main', rsid) {

  #Plotting
  p1 <- ggplot(whole_locus, aes(x = as.numeric(BP), y = -log10(P))) +
    geom_point(color = point_color) +
    geom_text_repel(data = whole_locus[c(whole_locus$P == min(whole_locus$P) | whole_locus$rsID %in% rsid), ],
              aes(x = as.numeric(BP), y = -log10(P), label = rsID),min.segment.length = 0) +
    theme_bw() +
    ylab(paste0("-log10(Pvalue)")) +
    xlab("") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())


  #Genes
  #Snp reference
  # we will need an additional mart for genes
  if (!exists("gene_ensembl")) {
    gene_ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = build)
  }
  #extract genes
  out_bm_genes_region <- getBM(
    attributes = c("start_position", "end_position", "ensembl_gene_id", "external_gene_name", "gene_biotype"),
    filters = c("chromosome_name", "start", "end"),
    values = list(str_replace(chr, "chr", ""), min, max),
    mart = gene_ensembl)
  ## define plot range for x-axis
  plot.range <- c(min(min, out_bm_genes_region$start_position),
                  max(max, out_bm_genes_region$end_position))
  ## rank gene_biotype label
  out_bm_genes_region_2 <- out_bm_genes_region
  tryCatch({
    out_bm_genes_region_2 <- out_bm_genes_region_2 %>%
      mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype), "protein_coding"),
             external_gene_name = fct_reorder2(external_gene_name, start_position, gene_biotype_fac, .desc = TRUE)) %>%
      filter(gene_biotype == "protein_coding")
  }, error = function(e) {
    message("Error: ", e)
  })

  #Plot genes
  random <- sample(1:nrow(out_bm_genes_region_2), nrow(out_bm_genes_region_2))
  p3 <- ggplot(data = out_bm_genes_region_2) +
    #geom_linerange(aes(x = random, ymin = start_position, ymax = end_position)) +
    coord_flip() +
    ylab(chr) +
    ylim(plot.range) +
    geom_text_repel(aes(x = 0, y = start_position, label = external_gene_name),
                   direction = 'y',max.overlaps = 100,min.segment.length = 100,
                   #segment.alpha=.2,segment.curvature=-.5,
              fontface = 2, alpha = I(0.7), hjust = "left", size = 2.5) +
    scale_y_continuous(n.breaks = 3) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.position = "bottom",
          panel.grid.major.y = element_blank())

  coloc.theme <- theme(plot.margin = margin(0, 0, 0, 0), panel.border = element_blank(),
                       panel.grid = element_blank(), axis.line = element_line(size = .2))
  p1 <- p1 + xlim(plot.range) + coloc.theme
  p3 <- p3 + coloc.theme + theme(axis.line.y  = element_blank())
  plot <- p1 / p3 + plot_layout(heights = c(3, 3))
  plot
}

#Plot
# we will need an additional mart for genes
gene_ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 38)
locuszoom(rsid = c(variant_split[,5],'rs34436714'))
ggsave(paste0(outdir,'figure2/nlrp12_locuszoom_minimal.pdf'),width = 3,height = 4)
