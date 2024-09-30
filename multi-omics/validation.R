suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(ComplexHeatmap)
  library(ggpubr)
})

covs <- "2000HIV/Phenotype/Phenotype_2000HIV_all_01.tsv"
outdir <- "2000HIV/cQTL/mofa/out/"
append <- "_corrected_scaled"
dir.create(paste0(outdir, "validation"))

#Load
model <- readRDS(paste0(outdir,"model", append, ".rds"))

#The weights
weights <- get_weights(model)

#Apply an extra regularization to the weights?
#Like, set all weights whose values are not significant? Would reduce overfitting and ease interpretation.
#X metabolites, Y proteins and Z methylation sites can reveal a factor related to A trait

alpha <- 0.01
#Add the mean and the sd of the weigths
weigths_df <- get_weights(model, as.data.frame = T) 
mean_w <- mean(weigths_df$value)
sd_w <- sd(weigths_df$value)
#Limit according to alpha
thr <- qnorm(alpha / 2, lower.tail = FALSE,
             mean = mean_w, sd = sd_w)

weights <- lapply(weights, function(x){
  y <- x
  #y <- apply(x, 2, scale) #Scaling by factor
  x[abs(y) < thr] <- 0
  return(x)
})
saveRDS(weights, paste0(outdir,"validation/validation_weights_alpha",alpha,append,".rds"))

eur <- fread(covs) %>% filter(ETHNICITY.White ==1) %>% pull(Record.Id)

#Read data
read_tsv <- function(file, name = "", correct = FALSE, norm = NULL, samples = "discovery", covdf, scale_it = T, ancestry = "White"){
      table <- fread(file)  %>% column_to_rownames(colnames(.)[1] )%>% t() %>% as.data.frame()

      if(samples == "discovery"){
        table <- table[!grepl("^ETZ", rownames(table)),]
      }else if(samples == "validation"){
        table <- table[grepl("^ETZ", rownames(table)),]
      }

      if(ancestry == "White"){
        table <- table[rownames(table) %in% eur,]
      }

      if(!is.null(norm)) {
        rn <- rownames(table)
        table <- apply(table, 2, norm) %>% as.data.frame()
        rownames(table) <- rn
      }



      table <- table %>% rownames_to_column("sample") %>% mutate(sample = gsub("OLVG","OLV",sample))
      if(correct){
        #Only the institutes
        institutedf <- model.matrix(as.formula("~-1+Institute.Abbreviation"),data=covdf[,c("sample","Institute.Abbreviation")])
        institutedf <- data.frame("sample" = covdf$sample, institutedf)
        table.corrected <- inner_join(table,institutedf, by= "sample")%>%
          mutate_at(colnames(table)[-1], function(f){
            #residuals(lm(f ~ age + sex + cd4 + cd8 + smoking + lockdown + season_sin + season_cos + ethnicityAsian +
            #               ethnicityBlack + ethnicityHispanic + ethnicityWhite, ., na.action=na.exclude))
            #residuals(lm(f ~ AGE + SEX_BIRTH + SMOKING + season_sin + season_cos + ETHNICITY.Asian +
            #               ETHNICITY.Black + ETHNICITY.White, ., na.action=na.exclude))
            residuals(lm(f ~ Institute.AbbreviationEMC + Institute.AbbreviationETZ + Institute.AbbreviationOLV +
                           Institute.AbbreviationRUMC, ., na.action = na.exclude))
          })
        table <- table.corrected[,colnames(table)]
      }
      #table <- melt(table, variable.name="feature")
      #if(name != "")table$view <- name
      if(scale_it) {
        rn <- table$sample
        table <- table %>% column_to_rownames("sample") %>% apply(2, scale, center = F) %>% as.data.frame()
        rownames(table) <- rn
      }else{
        table <- column_to_rownames(table, "sample")
      }
      return(table)
}


factor_list <- parallel::mclapply(c("discovery","validation"), function(cohort){
 
  correct <- ifelse(cohort == "discovery", T, F)

  #Correct covariates
  covdf1 <- fread(covs) %>% dplyr::select(-c(DNAm_SamplePlate,WGS,Cohort,ETHNICITY.Pacific_Islander,ETHNICITY.Native_American))%>%
    dplyr::rename("sample" = Record.Id)

  #Read files
  cyt <- read_tsv("meta_cQTL/data/2000HIV/Phenotype/cytokines_modified.tsv",
                  "cyt", correct = correct, norm = log1p, samples = cohort, covdf = covdf1)
  #RPMIs out
  cyt <- cyt[,!grepl("rpmi",colnames(cyt))]%>% t()
  #cyt %>% group_by(feature) %>% summarise(x=sum(!is.na(value))) %>% view()

  metab <- read_tsv(paste0("meta_cQTL/data/2000HIV-EU-",cohort,"/metabolites/phenotype.tsv"),
                    "metab", correct = correct, norm = log1p, samples = cohort, covdf = covdf1)%>%
    t()
  prot <- read_tsv(paste0("meta_cQTL/data/2000HIV-EU-",cohort,"/proteins/phenotype.tsv"),
                   "prot", correct = correct, samples = cohort, covdf = covdf1)%>%
    t()
  gex <- read_tsv(paste0("meta_cQTL/data/2000HIV-EU-",cohort,"/expression/phenotype.tsv"),
                 "gex", correct = correct, norm = log1p, samples = cohort, covdf = covdf1)%>%
    t()

  meth_file <- paste0(outdir, "meth.csv")
  meth <- read_tsv(meth_file, "meth", correct = correct, samples = cohort, covdf = covdf1)%>%
    t()


  #All rows and cols should be the same. Not gex cause I changed the names
  all(rownames(cyt) == rownames(weights$cyt))
      
  #Only common cyts, back to old cyts
   rownames(cyt) <- gsub("tnfa", "tnf", rownames(cyt)) %>%
      gsub("candconidia", "calbcon", .) %>% 
      gsub("ifng", "ifny", .) %>%
      gsub("candhyphae", "calbhy", .) 
    cyt <- cyt[rownames(weights$cyt),]

  all(rownames(metab) == rownames(weights$metab))
  all(rownames(prot) == rownames(weights$prot))
  all(rownames(meth) == rownames(weights$meth))

      data_list_log <- list("cyt" = cyt, "metab" = metab, "prot" = prot, "gex" = gex, "meth" = meth)
      data_list <- lapply(data_list_log, function(x) {
        out <- apply(x,2,scale)
        rownames(out) <- rownames(x)
        return(out)
        })

      factors <- lapply(names(data_list), function(n){
        print(n)
        data_clean <- apply(data_list[[n]],1,function(x) ifelse(is.na(x), 0, x))
        res <- data_clean %*% weights[[n]]
      })

      #Common ppl
      factors_int <- lapply(factors, function(x) x[Reduce(intersect, lapply(factors, rownames)),])
      factors_sum <- Reduce("+", factors_int) %>% apply(2,scale)
      rownames(factors_sum) <- rownames(factors_int[[1]])

  #Save the input data
  if (cohort == "validation") {
    lapply(names(data_list_log), function(d){
      fwrite(as.data.frame(data_list_log[[d]][,rownames(factors_sum)]),
             paste0(outdir, "validation/",d,"_input",append,".csv"), row.names = T)
  })
    #Save the validation values
   fwrite(as.data.frame(factors_sum), paste0(outdir, "validation/validation_factors_alpha",alpha, append, ".csv"), row.names = T)
        #Save covariates
   fwrite(covdf1 %>% filter(sample %in% rownames(factors_sum)),
          paste0(outdir, "validation/validation_covs.csv"), row.names = T)
  }

  return(factors_sum)

}, mc.cores = 2)

saveRDS(factor_list,
        paste0(outdir, "validation/factor_list_alpha", alpha, append, ".rds"))

#Check discovery
#1. melt matrices to a "factor" column,
#merge by sample id to mofa factors, correlate

factor_disc <- factor_list[[1]] %>% as.data.frame() %>%
  rownames_to_column("id") %>%
  melt(id.vars = "id", variable.name = "factor", value.name = "calc")

mofa_factors <- get_factors(model)%>% as.data.frame() %>%
rownames_to_column("id")%>% melt(id.vars = "id", variable.name = "factor", value.name = "mofa")%>%
mutate("factor" = gsub("single_group.", "", factor))

mofa_and_calculated <- inner_join(factor_disc, mofa_factors)

ggplot(mofa_and_calculated, aes(x=calc, y=mofa))+
  geom_point() + facet_wrap(~factor)+
  stat_cor() + geom_smooth(method = "lm")+
  theme_bw() + theme(panel.grid = element_blank())

ggsave(paste0(outdir,
              "validation/mofa_and_calculated_cor_alpha", alpha, append, ".pdf"))
