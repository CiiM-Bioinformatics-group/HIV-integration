suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
    library(ggpubr)
    library(MOFA2)
    library(reticulate)
    library(clusterProfiler)
    library(parallel)
  INT <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
})

reticulate::use_condaenv('mofa')
reticulate::use_python('~/miniconda3/envs/mofa/bin/python')

covs = '2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv'
outdir = "2000HIV/cQTL/mofa/out/"
dir.create(outdir, showWarnings = F)
dir.create(paste0(outdir,'data'), showWarnings = F)

eur <- fread(covs) %>% filter(ETHNICITY.White ==1) %>% pull(Record.Id)

mclapply(c(T,F), function(correct){
  #Reading the matrixqtl files
  read_tsv <- function(file, name = '', correct = FALSE, norm = NULL, samples = 'discovery', ancestry = 'White'){
    table <- fread(file)  %>% column_to_rownames(colnames(.)[1] )%>% t() %>% as.data.frame()

    if(samples == 'discovery'){
      table <- table[!grepl('^ETZ', rownames(table)),]
    }else if(samples == 'validation'){
      table <- table[grepl('^ETZ', rownames(table)),]
    }
    if(ancestry == 'White'){
      table <- table[rownames(table) %in% eur,]
    }
    if(!is.null(norm)) table <- apply(table, 2, norm) %>% as.data.frame()

    table <- table %>% rownames_to_column('sample') %>% mutate(sample = gsub('OLVG','OLV',sample))
    if(correct){
      #Only the institutes
      institutedf <- model.matrix(as.formula("~-1+Institute.Abbreviation"),data=covdf[,c('sample','Institute.Abbreviation')])
      institutedf <- data.frame('sample' = covdf$sample, institutedf) %>%
        inner_join(covdf[,c('sample','AGE','SEX_BIRTH')], by = 'sample')
      table.corrected <- inner_join(table,institutedf, by= 'sample')%>%
        mutate_at(colnames(table)[-1], function(f){
          residuals(lm(f ~ AGE + SEX_BIRTH + Institute.AbbreviationEMC + Institute.AbbreviationETZ + Institute.AbbreviationOLV +
                          Institute.AbbreviationRUMC, ., na.action = na.exclude))
        })
      table <- table.corrected[,colnames(table)]
    }
    table <- melt(table, variable.name='feature')
    if(name != '')table$view <- name
    return(table)
  }

  #Correct covariates
  covdf <- fread(covs) %>% dplyr::select(-c(DNAm_SamplePlate,WGS,Cohort,ETHNICITY.Pacific_Islander,ETHNICITY.Native_American))%>%
    dplyr::rename('sample' = Record.Id)

  #Read files
  cyt <- read_tsv('meta_cQTL/data/2000HIV/Phenotype/cytokines_modified.tsv',
                  'cyt', correct = correct, norm = log1p)
  #RPMIs out
  cyt <- cyt %>% subset(!grepl('rpmi',feature))
  #cyt %>% group_by(feature) %>% summarise(x=sum(!is.na(value))) %>% view()

  metab <- read_tsv('meta_cQTL/data/2000HIV-EU-discovery/metabolites/phenotype.tsv',
                    'metab', correct = correct, norm = log1p)
  prot <- read_tsv('meta_cQTL/data/2000HIV-EU-discovery/proteins/phenotype.tsv',
                   'prot', correct = correct)
  gex <- read_tsv('meta_cQTL/data/2000HIV-EU-discovery/expression/phenotype.tsv',
                  'gex', correct = correct, norm = log1p)

  meth_file <- paste0(outdir, 'meth.csv')
  meth <- read_tsv(meth_file, 'meth', correct = correct)

  #data_list <- list('cyt' = cyt, 'metab' = metab, 'prot' = prot, 'gex' = gex, 'meth' = meth)
  #lapply(names(data_list), function(d){
  #  samples <- Reduce(intersect, sapply(data_list, colnames))
  #  fwrite(as.data.frame(data_list_log[[d]][,rownames(factors_sum)]), paste0(outdir, 'data/',d,'_input',
  #                                                                       ifelse(correct,'_corrected',''),'.csv'), row.names = T)
  #})

  #Save all

  all <- rbind(cyt,metab,prot,gex,meth)

  #Only samples with all modalities
  full_samples <- Reduce(intersect,list(unique(cyt$sample), unique(metab$sample),
                            unique(prot$sample),unique(gex$sample),unique(meth$sample)))

  all <- filter(all, sample %in% full_samples)

  #MOFA
  MOFAobject <- create_mofa(all)
  plot_data_overview(MOFAobject)


  mclapply(c(T,F), function(scale.it){
    #Options
    data_opts <- get_default_data_options(MOFAobject)
    model_opts <- get_default_model_options(MOFAobject)
    train_opts <- get_default_training_options(MOFAobject)
    model_opts$num_factors <- 30
    train_opts$verbose <- TRUE
    data_opts$scale_views <- scale.it
    #Prepare model
    MOFAobject <- prepare_mofa(
      object = MOFAobject,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )

    if(correct & scale.it){
       outfile = paste0(outdir,"model_corrected_scaled.hdf5")
     }else if(correct){
       outfile = paste0(outdir,"model_corrected_unscaled.hdf5")
     }else if(scale.it){
       outfile = paste0(outdir,"model_uncorrected_scaled.hdf5")
     }else{
       outfile = paste0(outdir,"model_uncorrected_unscaled.hdf5")
     }
    train_opts$convergence_mode < 'slow'
    MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = FALSE)
  }, mc.cores = 2)

}, mc.cores = 4)
