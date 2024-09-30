#Writing a table with the selected clinical variables and the factors calculated
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(MOFA2)
  library(parallel)
})

outdir = "2000HIV/cQTL/mofa/out/"
dir.create(paste0(outdir,'clin_cor'))

#mclapply(c('corrected','uncorrected'), function(i){
#  mclapply(c('scaled','unscaled'), function(j){
#    append = paste0('_',i,'_',j)
append = '_corrected_scaled'
    #Load
    model <- readRDS(paste0(outdir,'model', append, '.rds'))
    factors <- do.call("rbind",get_factors(model))

    #Clinical data
    clin <- '2000HIV/Phenotype/220317_2000hiv_study_export_processed.csv'
    clin2 <- '2000HIV/Phenotype/231912_CVD_phenotypes.csv'
    info <- '2000HIV/Phenotype/220221_ECRF_GUIDE_2000HIV.docx'

    clin_df <- fread(clin)%>%as.data.frame()%>%inner_join(fread(clin2)%>%dplyr::rename('Record Id' = `Record.Id`)%>%
                                                            .[,-c(2:6,8,9)], by = 'Record Id')%>%
      inner_join(as.data.frame(factors) %>% rownames_to_column('Record Id'))
    
    colnames(clin_df) <- gsub('#', '%', colnames(clin_df))

    #Manually select the variables
    #All the MED%: "Use of medication the previous two weeks besides antiretroviral medication (ART)"
    #All the CVDRISC%: "Cardiovascular risk factors"
    #All the MH_SOA%: "Did patient ever have one of these sexual transmittable infections? (STIs/STDs)"
    #All the OPP_PRECART%: "Where there any opportunistic infections before start ART?"
    #All the CART%: "CART	ART regime at the time of baseline visit"
    #All the CART_NRTI%: "Specify NRTI"
    #All the CART_NTRTI%: "Specify NTRTI"
    #All the CART_NNRTI%: "Specify NNRTI"
    #All the CART_II%: "Specify Integrase inhibitors"
    #HIV_STAGE_CDC: "CDC stage at the time of start ART"
    #VL_LATEST: "What is the most recent viral load value at baseline?"
    #All the OPP_POSTCART%: "Where there any opportunistic infections after start ART?"
    #MALIG POSTCART: "Where there any malignancy after start ART?"
    #RP_YESNO: "Rapid progressors" #SET 2 TO MISSING! #TODO #MODIFY
    #INR_YESNO: "Immunological non responders" #SET 2 TO MISSING! #TODO #MODIFY
    #CONTROLLER: "Elite controller" #2 means elite here, not missing
    #ACUTE_HIV: "Everything that is not 0 means the HIV infection is recent" #SET 4 TO MISSING #TODO #MODIFY
    #MH_TR_END and MH_TR_END_D%: "Endocrine (hormonal) or metabolic disorders diagnosed"
    #MH_TR_CV and MH_TR_CV_D%: "Cardiovascular or thromboembolic disorders diagnosed?"
    #MH_TR_RESP and MH_TR_RESP_D%: "Respiratory disorders  diagnosed?"
    #MH_TR_GE and MH_TR_GE_D%: "Gastrointestinal disorders diagnosed?"
    #HEPA_BASELINE: "Has patient ever had hepatitis A (HAV) disease or a vaccination?": #TODO Change to factor
    #1 is previous, 2 is vaccine, 3 is unknown, 5 is ig present
    #HEPB_BASELINE: "Has patient ever had hepatitis B (HBV) disease or a vaccination?": #TODO Change to factor
    #1 is previous, 2 is vaccine, 3 is unknown, 5 is ig present
    #HEPC: "HCV in medical history"
    #TOXO%: "Has patient ever been tested for toxoplasm and what was the result?"
    #CMV%: "Has patient ever been tested for CMV and what was the result?"
    #Select those and make the changes

    grep_col <- function(df,col) {
      if(!any(grepl(col,colnames(df)))) stop(paste(col, 'is not in the colnames'))
      colnames(df)[grepl(col,colnames(df))]
    }
    clin_df_clean <- clin_df %>%
      column_to_rownames('Record Id')%>%
      dplyr::select(
        #grep_col(.,'^MED%'), grep_col(.,'^CVDRISC%'),
        #grep_col(.,'^MH_SOA%'), grep_col(.,'^OPP_PRECART%'), grep_col(.,'^CART%'), grep_col(.,'^CART_NRTI%'),
        #grep_col(.,'^CART_NTRTI%'), grep_col(.,'^CART_NNRTI%'), grep_col(.,'^CART_II%'), grep_col(.,'^HIV_STAGE_CDC$'),
        #grep_col(.,'^VL_LATEST$'), grep_col(.,'^OPP_POSTCART%'), grep_col(.,'^MALIG_POSTCART$'),
        #grep_col(.,'^VL_LATEST_OTHER$'),grep_col(.,'^IMT_$'),
        grep_col(.,'^RP_YESNO$'), grep_col(.,'^INR_YESNO$'),
        #grep_col(.,'^CONTROLLER$'), grep_col(.,'^ACUTE_HIV$'), grep_col(.,'^HEPA_BASELINE$'), 
        #grep_col(.,'^CVD2y$'),grep_col(.,'^CVDdenovo$'),
        grep_col(.,'PLAQUE'),
        grep_col(.,'^MH_TR_END$'), grep_col(.,'^MH_TR_END_D%'),  grep_col(.,'^MH_TR_CV$'), grep_col(.,'^MH_TR_CV_D%'),
        grep_col(.,'^MH_TR_RESP$'), grep_col(.,'^MH_TR_RESP_D%'), grep_col(.,'^MH_TR_GE$'), grep_col(.,'^MH_TR_GE_D%'),
        #grep_col(.,'^MH_TR_PSY$'), grep_col(.,'^MH_PSY_D%'), grep_col(.,'^MH_TR_LOCO$'), grep_col(.,'^MH_TR_LOCO_D%'),
        #grep_col(.,'^MH_SURG$'), grep_col(.,'^MH_SURG_D%'),
        #grep_col(.,'^HEPB_BASELINE$'), grep_col(.,'^HEPC$'),
        #grep_col(.,'^TOXO%'), grep_col(.,'^CMV%'),
        grep_col(.,'^Factor')
      )%>%
      #Necessary changes
      mutate(RP_YESNO = ifelse(RP_YESNO == 2, NA, RP_YESNO),
             INR_YESNO = ifelse(INR_YESNO == 2, NA, INR_YESNO),
             #ACUTE_HIV = ifelse(ACUTE_HIV == 4, NA, ACUTE_HIV),
             #ACUTE_HIV = ifelse(ACUTE_HIV != 0, 1, ACUTE_HIV),
             #HEPA_BASELINE = as.character(HEPA_BASELINE),
             #HEPB_BASELINE = as.character(HEPB_BASELINE),
             )

    fwrite(clin_df_clean, paste0(outdir,'clin_cor/clin_df',append,'.csv'), row.names=T)

#  }, mc.cores=2)
#}, mc.cores=4)
