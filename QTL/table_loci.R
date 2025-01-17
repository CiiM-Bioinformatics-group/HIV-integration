suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(ieugwasr)
# library(phenoscanner)
  library(magrittr)
  library(plinkbinr)
})

# palette from https://github.com/jbgb13/peRReo/blob/master/R/latinpalettes.R
palette <- c("#7E0734","#04AAB3","#0F3645","#603A79","#A0405D","#246B82","#B3B4AF","#F5F561","#4B1E3D")

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
clump_kb <- args[2]
clump_r2 <- args[3]
outfile <- args[4]
plotdir <- args[5]
cov <- args[6]
mapping <- args[7]
bfile <- args[8]
num_samples_file <- args[9]

#Get prefix bfile
bfile <- gsub(".bed$", "", bfile)

message("input ", input,
        "\nclump_kb ", clump_kb,
        "\nclump_r2 ", clump_r2,
        "\noutfile ", outfile,
        "\nplotdir ", plotdir,
        "\ncov ", cov,
        "\nmapping ", mapping,
        "\nbfile ", bfile,
        "\nnum_samples_file ", num_samples_file)

stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

plotdir <- paste0(plotdir, "/loci")
if (!dir.exists(plotdir)) dir.create(plotdir)

#Load genomewide
message("Loading genome-wide SNPs...")
gw <- fread(input)

# extra pval filter if pval exists
if (exists("pval")) {
  gw <- gw %>% 
    mutate(`p-value` = as.numeric(`p-value`)) %>% 
    filter(`p-value` < pval)
}

if (nrow(gw) == 0) {
  message("No genomewide significant snps found")
  fwrite(as.data.frame(gw), outfile, sep = "\t")
} else {
  loci <- gw %>%
    dplyr::rename("pval" = `p-value`, "rsid" = SNP) %>%
    ld_clump_local(clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = 1, 
      bfile = bfile, plink_bin = plinkbinr::get_plink_exe()) %>%
    #mutate(CHR = str_replace(CHR, "^chr|^CHR", "")) %>%
    dplyr::rename("SNP" = rsid) %>%
    tidyr::separate(SNP, into = c("SNPid", "rsid"), sep = ";") %>%
    tidyr::separate(SNPid, into = c("CHR", "BP", "REF", "ALT"), remove = F) %>%
    arrange(factor(CHR, levels = paste0("chr", 1:22)), as.numeric(BP)) %>%
    rownames_to_column("locus") %>% # this is just temporary to get it in first pos
    setDT() 

  # give each individual snp (locus) a unique id
  loci[, "locus" := .GRP, by = SNPid]

  loci <- loci %>% dplyr::select(-SNPid)

  if (FALSE) { # phenoscanner dead
    #Add phenoscanner info about gwas, eqtl and annotation
    ps <- phenoscanner(loci$rsid, catalogue = "GWAS")
    if (nrow(ps$snps) > 0) {
      pos <- ps$snps %>%
        select(c("snp", "consequence", "hgnc")) %>%
        rename("rsid" = snp)
      loci %<>%
        left_join(pos, "rsid")
    }
    if (nrow(ps$results) > 0) {
      gwas <- ps$results %>%
        select(c("snp", "trait")) %>%
        rename("rsid" = snp)
      loci %<>%
        left_join(gwas, "rsid") %>%
        group_by(across(c(-trait))) %>% # at this point t-stat^1 and consequence??2 are created
        summarise_at("trait", function(x) paste(unique(x), collapse = ";"))
    }

    eqtl <- phenoscanner(loci$rsid, catalogue = "eQTL")$results
    if (nrow(eqtl) > 0) {
      eqtl <- eqtl %>%
        dplyr::select(c("rsid", "exp_gene"))

      loci %<>%
        left_join(eqtl, "rsid") %>%
        group_by(across(c(-exp_gene))) %>%
        summarise_at("exp_gene", function(x) paste(unique(x), collapse = ";")) %>%
        mutate("exp_gene" = gsub("-", "", exp_gene))
    }
  }

  loci %<>%
    arrange(as.numeric(locus)) %>%
    ungroup() #%>%
    #dplyr::select(-id)

  #loci_cp <- loci

  # Add the interactions if we are creating a table of the main effects.
  #if (cov == "main") {
  if (FALSE) {
    tryCatch({
      loci_interaction <- lapply(seq_len(nrow(loci)), function(row) {
        dat <- loci[row, ]
        #Get all mapping results for the gene
        subdirs <- list.dirs(paste0(mapping, "/", dat$CHR))
        files <- lapply(subdirs, function(subdir) {
          file_path <- file.path(subdir, paste0(dat$gene, ".txt"))
          if (file.exists(file_path) && !endsWith(subdir, "/main")) return(file_path)
        })
        files <- Filter(Negate(is.null), files)  # remove the returned "NULL" elements

        # Check if there actually files
        if (length(files) == 0) stop()

        # Get p-values
        df <- lapply(files, function(f) {
          interact <- gsub(mapping, "", f) %>%
            gsub(dat$gene, "", .) %>%
            gsub(dat$CHR, "", .) %>%
            gsub(".txt", "", .) %>%
            gsub("/", "", .)
          df <- fread(cmd = paste0("grep -w '", dat$CHR, ":", dat$BP, ":",
            dat$REF, ":", dat$ALT, ";", dat$rsid, "' ", f))
          rdf <- data.frame(df$V5)
          colnames(rdf) <- interact
          return(rdf)
        }) %>%
        bind_cols()

        df$rsid <- dat$rsid
        df$gene <- dat$gene
        return(df)
      }) %>%
      bind_rows()

      loci <- loci %>%
        left_join(loci_interaction, by = c("rsid", "gene"))
    }, error = function(e) {
      print(e)
      print("No interactions found")
    })
  }

  # add num of samples for each phenotype
  num_samples <- fread(num_samples_file)
  # if just one line, add the number of samples to all loci
  if (nrow(num_samples) == 1) {
    loci$num_samples <- num_samples$num_samples[1]
  } else {
    loci <- loci %>%
      left_join(num_samples, by = c("gene" = pheno))
  }

  fwrite(as.data.frame(loci), outfile, sep = "\t")
  message(outfile)

  if (FALSE) {
    tryCatch({
      if (all(grepl("^pbmc_|^wb_", loci$gene))) {
        #Plot distribution of some variables only if the name corresponds to the standard
        lapply(c("cytokine", "stimulation", "consequence"), function(y) {
          loci %>%
            tidyr::separate(gene, into = c("cell", "time", "cytokine", "stimulation"), remove = FALSE) %>%
            mutate(gene = toupper(gene) %>%
            gsub("IL", "IL-", .) %>%
            gsub("G", "g", .) %>%
            gsub("A", "a", .)) %>%
            dplyr::rename(yaxis = y) %>%
            group_by(yaxis) %>%
            summarise(n = n()) %>%
            ggplot(aes(x = n, y = reorder(yaxis, n), fill = yaxis)) +
            geom_col() +
            geom_text(aes(label = n), color = "white", hjust = 0, nudge_x = -.5, size = 4, fontface = "bold") +
            scale_fill_manual(values = grDevices::colorRampPalette(palette[1:9])(10)) +
            coord_cartesian(clip = "off") +
            scale_x_continuous(expand = c(.01, .01)) +
            theme_void() +
            theme(
              axis.text.y = element_text(size = 14, hjust = 1),
              ## make sure labels doesn"t get cut, part 2
              plot.margin = margin(15, 30, 15, 15),
              legend.position = 0, strip.placement = "outside"
            ) +
            ggtitle(paste0("Independent loci per\n", y))
          plotout <- paste0(plotdir, "/", cov, "_", y, ".pdf")
          message(plotout)
          ggsave(plotout, width = 4, height = 4)
        })
      }
    }, error = function(e) {
      message("plotting failed ", e)
    })
  }
}




