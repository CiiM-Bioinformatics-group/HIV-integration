# Readme
# Created on: 26 October 2022
# Source of raw genotype data: /mnt/workspace/1. Original data/1.Raw/2000HIV/2000HIV_GWAS_GSAData/GSA796/PLINK_220222_0252
# File names: GSA796.phenotype, GSA796.ped, GSA796.map, GSA796.bat, GSA796.script

# Part of analysis (QC) was done on DRE: data/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/
# PCA analysis and alignement to reference genome was done locally

# Step 1 - Convert to bed files ------------------------------------------------------------------------------------------------------------------------------------------
# Options in effect:
  --file GSA796
  --make-bed
  --out 2000HIV_GSA796_all_chromosomes


# Step 2 - Find the duplicated samples (OLV282,OLV283,OLV284) and convert to unique ids ------------------------------------------------------------------------------------------------------------------------------------------
# The duplicated values (IID) change from "1134" "982"  "802" to "1134.1" "982.1"  "802.1"  in the fam file
# Check R script: 1_find_duplicated_samples_in_fam_file_29nov2022.R

# Step 3 - Calculate missingness per sample and SNP using the updated fam file (2000HIV_GSA796_all_chromosomes.fam) from step 2 ------------------------------------------------------------------------------------------------------------------------------------------
# Options in effect:
--bfile 2000HIV_GSA796_all_chromosomes
--missing
--out miss_2000HIV_GSA796_all_chromosomes

# output
# sample missing data report written to 
# miss_2000HIV_GSA796_all_chromosomes.imiss
# variant-based missing data report written to
# miss_2000HIV_GSA796_all_chromosomes.lmiss

# Investigate missingness per individual and per SNP and make histograms.
# See R script: 2000HIV_GWAS_QC_plots_29nov2022.R

# Find duplicated samples with high missingness
# Check R script: 1_find_duplicated_samples_in_fam_file_sex_and_center_info_29nov2022.R

# Step 4 - Remove duplicated samples ------------------------------------------------------------------------------------------------------------------------------------------
# Options in effect:
--bfile 2000HIV_GSA796_all_chromosomes
--make-bed 
--remove samples_with_higher_F_MISS_to_remove_29nov2022.txt
--out new_2000HIV_GSA796_all_chromosomes

# Step 5 - Keep autosomes and sex chromosomes
# The autosomes should be coded 1 through 22. The following other codes can be used to specify other chromosome types:
#     X    X chromosome                    -> 23
#     Y    Y chromosome                    -> 24
#     XY   Pseudo-autosomal region of X    -> 25
#     MT   Mitochondrial                   -> 26

# Options in effect:
--bifle new_2000HIV_GSA796_all_chromosomes
--chr 1-24
--make-bed --out new2_2000HIV_GSA796_all_chromosomes

# Step 5 - Update sex in the plink files ------------------------------------------------------------------------------------------------------------------------------------------
# Note
# On the study export (castor) 0=male and 1= female (DRE path: data/1. Original data/2.Processed/2000HIV/2000HIV_ClinicalData_eCRFs_CASTOR/221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.csv)
# Plink files use: 1/M = male, 2/F = female, 0 = missing
# Therefore, we coded females with 2 (n 280), males with 1 (n 1573) and NAs (n = 61) with 0

# Options in effect:
--bifle new2_2000HIV_GSA796_all_chromosomes
--make-bed
--update-sex to_update_sex_for_new_2000HIV_GSA796_all_chromosomes_29nov2022.txt
--out updated_sex_new2_2000HIV_GSA796_all_chromosomes

# Step 6 - Check for sex discrepancies ------------------------------------------------------------------------------------------------------------------------------------------

# Options in effect:
--bfile updated_sex_new2_2000HIV_GSA796_all_chromosomes
--check-sex
--out check_sex_updated_sex_new2_2000HIV_GSA796_all_chromosomes
# 722909 variants and 1914 people


# QC per marker --------------------------------------------------------------------------------------------------------------------------------------
# Step 7 - Remove SNPs with a call rate < 95%
# All ethnicities
# Options in effect 
--bfile updated_sex_new2_2000HIV_GSA796_all_chromosomes 
--geno 0.05 # 2720 variants removed due to missing genotype data
--make-bed 
--out all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796

# 720189 remained

# Step 8 - Remove SNPs that fail HWE < 1e-6 ------------------------------------------------------------------------------------------------------------------------------------------
# Testing genetic markers for Hardy–Weinberg equilibrium (HWE) is an important tool for detecting genotyping errors in large-scale genotyping studies. 
# For markers at the X chromosome, typically the χ2 or exact test is applied to the females only, and the hemizygous males are considered to be uninformative.

# All ethnicities  
# Exclude X chromosome 
# Options in effect 
--bfile all_ethnicities_geno0.05_european_updated_sex_new2_2000HIV_GSA796 
--chr 1-22
--make-bed 
--out all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes # 686708 varians on autosomes

# Remove autosomal SNPs that fail HWE at P < 1e-6 using only samples of white ancestry
# Options in effect 
--bfile all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
--hwe 1e-6 # 471 variants removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_whites.txt
--out whites_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes

# Remove autosomal SNPs that fail HWE at P < 1e-6 using only samples of black ancestry
# Options in effect 
--bfile all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
--hwe 1e-6 # 50 variants removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_blacks.txt
--out blacks_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes

# Remove autosomal SNPs that fail HWE at P < 1e-6 using only samples of Asian ancestry
# Options in effect 
--bfile all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
--hwe 1e-6 # 11 variants removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_asians.txt
--out asians_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes

# Remove autosomal SNPs that fail HWE at P < 1e-6 using only samples of Hispanic ancestry
# Options in effect 
--bfile all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
--hwe 1e-6 # 2 variants removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_hispanics.txt
--out hispanics_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes

# Remove autosomal SNPs that fail HWE at P < 1e-6 using only samples of Native American ancestry
# Options in effect 
--bfile all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
--hwe 1e-6 # 0 variants removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_native_american.txt
--out native_americans_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes

# Remove autosomal SNPs that fail HWE at P < 1e-6 using only samples of mixed ancestry
# Options in effect 
--bfile all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
--hwe 1e-6 # 19 variants removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_mixed.txt
--out mixed_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes


# Keep only X chromosome and select females (n = 28)
# Options in effect 
--bfile all_ethnicities_geno0.05_european_updated_sex_new2_2000HIV_GSA796 
--chr 23 # 28526 variants on X chromosome
--keep only_females_for_all_ethnicities_cohort_1dec2022.txt # 307 European females # This file was creaed using script: 2_2000HIV_GWAS_QC_plots_29nov2022.R
--make-bed 
--out all_ethnicities_females_geno0.05_updated_sex_new2_2000HIV_GSA796_X_chromosome

# Remove SNPs that fail HWE at P < 1e-6 on X chromosome in samples of white ancestry
# Options in effect 
--bfile all_ethnicities_females_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome
--hwe 1e-6 #  0 variants on X chromosome removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_whites.txt
--out whites_females_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome

# Remove SNPs that fail HWE at P < 1e-6 on X chromosome in samples of Black ancestry
# Options in effect 
--bfile all_ethnicities_females_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome
--hwe 1e-6 #  1 variants on X chromosome removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_blacks.txt
--out blacks_females_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome

# Remove SNPs that fail HWE at P < 1e-6 on X chromosome in samples of Asian ancestry
# Options in effect 
--bfile all_ethnicities_females_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome
--hwe 1e-6 #  0 variants on X chromosome removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_asian.txt
--out asians_females_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome

# Remove SNPs that fail HWE at P < 1e-6 on X chromosome in samples of Hispanic ancestry
# Options in effect 
--bfile all_ethnicities_females_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome
--hwe 1e-6 #  0 variants on X chromosome removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_hispanic.txt
--out hispanics_females_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome

# Remove SNPs that fail HWE at P < 1e-6 on X chromosome in samples of Native American ancestry
# Options in effect 
--bfile all_ethnicities_females_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome
--hwe 1e-6  
--make-bed
--keep 2000HIV_native_american.txt
--out native_americans_females_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome

# No people remaining after --keep. No females of Native American ancestry exist

# Remove SNPs that fail HWE at P < 1e-6 on X chromosome in samples of Mixed ancestry
# Options in effect 
--bfile all_ethnicities_females_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome
--hwe 1e-6 #  0 variants on X chromosome removed due to Hardy-Weinberg exact test
--make-bed
--keep 2000HIV_mixed.txt
--out mixed_females_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome


# Exclude SNPs from X chromosome and autosomes that failed HWE
# Options in effect 
--bfile all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796 
--exclude snps_that_failed_hwe_in_all_ethnicities_in_autosomes_and_X_chrom_2dec2022.txt # This file was created using script: 2_2000HIV_GWAS_QC_plots_29nov2022.R
--chr 1-23
--make-bed
--out all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
# 714696 variants and 1914 samples (all ethnicities)
# Warning: 26301 het. haploid genotypes present (see: discovery_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_chr1to23.hh)


# Step 9 ------------------------------------------------------------------------------------------------------------------------------------------
# Given that different ancestries have different MAF, we tested to apply a MAF filter per ethnicity group

# Remove SNPs with MAF < 1% in samples of white ancestry
# Options in effect 
--bfile all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--maf 0.01 # 517432 variants remain; 197264 SNPs were removed due to minor allele thresholds
--make-bed 
--keep 2000HIV_whites.txt # 1381 samples of white ancestry remain
--out whites_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 

# Remove SNPs with MAF < 1% in samples of black ancestry
# Options in effect 
--bfile all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--maf 0.01 # 411066 variants remain; 303630 SNPs were removed due to minor allele thresholds
--make-bed 
--keep 2000HIV_blacks.txt # 186 samples of black ancestry remain
--out blacks_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 

# Remove SNPs with MAF < 1% in samples of Asian ancestry
# Options in effect 
--bfile all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--maf 0.01 # 431940 variants remain; 282756 SNPs were removed due to minor allele thresholds
--make-bed 
--keep 2000HIV_asians.txt # 85 samples of black ancestry remain
--out asians_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 

# Remove SNPs with MAF < 1% in samples of Hispanic ancestry
# Options in effect 
--bfile all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--maf 0.01 # 503718 variants remain; 210978 SNPs were removed due to minor allele thresholds
--make-bed 
--keep 2000HIV_hispanics.txt # 48 samples of hispanic ancestry remain
--out hispanics_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 

# Remove SNPs with MAF < 1% in samples of Native american ancestry
# Options in effect 
--bfile all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--maf 0.01 # 255236 variants remain; 459460 SNPs were removed due to minor allele thresholds
--make-bed 
--keep 2000HIV_native_americans.txt # 3 samples of Native American ancestry remain
--out native_americans_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 

# Remove SNPs with MAF < 1% in samples of Mixed ancestry
# Options in effect 
--bfile all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--maf 0.01 # 489194 variants remain; 225502 SNPs were removed due to minor allele thresholds
--make-bed 
--keep 2000HIV_mixed.txt # 128 samples of black ancestry remain
--out mixed_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 

# Step 11 - QC per sample ------------------------------------------------------------------------------------------------------------------------------------------

# Exclude samples with a call rate < 0.975 (97.5%) in the cohort of all ethnicities
# Options in effect 
--bfile all_ethnicities_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--mind 0.025
--make-bed
--out all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
#  714696 variants and  1903 people pass filters,  11 samples were removed due to missing genotype data
# total genotyping rate in remaining samples is 0.998805

# Step 12 - Remove individuals with a heterozygosity rate deviating more than 3 sd from the mean --------------------------------------------------------------------------------------------------------------------------------------------------------------
# Identify individuals with excess heterozygosity rate per ethnicity group
# Checks for heterozygosity are performed on a set of SNPs which are not highly correlated.
# Therefore, to generate a list of non-(highly)correlated SNPs, we exclude high LD regions () and prune the SNPs using the command --indep-pairwiseí.
# The parameters 50 5 0.2 stand respectively for: the window size, the number of SNPs to shift the window at each step, and the multiple correlation coefficient for a SNP being regressed on all other SNPs simultaneously.
# Source of high LD regions: https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt

# Important note - Exclude the X chromosome to calculate heterozygosity rate

# Heterozygosity Excess  is an indicator of the quantity of excess heterozygote calls relative to expectations based on Hardy Weinberg Equilibrium (HWE).
# It ranges from -1 (complete deficiency of heterozygotes) to 1 (100% heterozygotes). Autosomal variants with values less than -0.4 or greater than 0.4 are expected to be of bad quality.
# Heterozygosity excess on the X-chromosome is dependent on the proportion of males in your study. For example, if your study consists of only males, X-chromosomal variants will 
# have a heterozygosity excess of -1.  

# Subset samples of white ancestry
# Options in effect 
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--keep 2000HIV_whites
--make-bed 
--chr 1-22
--exclude high-LD-regions-hg19-GRCh37.txt # 20777 variants excluded and 665394 variants remaining
--range 
--indep-pairwise 50 5 0.2 
--out whites_indepSNP
# 412558 of 665394 variants removed

# Extract the pruned SNPs
# Options in effect 
--bfile  all_ethnicities_mind0.025_maf0.01_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--extract whites_indepSNP.prune.in
--keep 2000HIV_whites
--chr 1-22
--make-bed 
--out whites_only_autosomes_pruned
# 252836 variants and 1373 people pass filters and QC

# Calculation of heterozygosity rate to identify individuals with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes
# Options in effect 
--bfile whites_only_autosomes_pruned
--het
--out whites_het_only_autosomes # Use this file to calculate heterozygosity rate
# Report written to whites_het_only_autosomes.het
# 32 samples of White ancestry were identified with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes

# Calculate missingness using only autosomes
# For white ancestry
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--chr 1-22
--missing
--keep 2000HIV_whites.txt
--out whites_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
# 686171 out of 714696 variants loaded from bim file

# Subset samples of Black ancestry
# Options in effect 
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--keep 2000HIV_blacks.txt
--make-bed 
--chr 1-22
--exclude high-LD-regions-hg19-GRCh37.txt # 20777 variants excluded and 665394 variants remaining
--range 
--indep-pairwise 50 5 0.2 
--out blacks_indepSNP
# 380212 of 665394 variants removed

# Extract the pruned SNPs
# Options in effect 
--bfile  all_ethnicities_mind0.025_maf0.01_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--extract blacks_indepSNP.prune.in
--keep 2000HIV_blacks.txt
--chr 1-22
--make-bed 
--out blacks_only_autosomes_pruned
# 285182 variants and 184 people pass filters and QC

# Calculation of heterozygosity rate to identify individuals with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes
# Options in effect 
--bfile blacks_only_autosomes_pruned
--het
--out blacks_het_only_autosomes # Use this file to calculate heterozygosity rate
# Report written to blacks_het_only_autosomes.het

#  samples of Black ancestry were identified with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes

# Calculate missingness using only autosomes
# For all ethnicities
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--chr 1-22
--missing
--keep 2000HIV_blacks.txt
--out blacks_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
# 686171out of 714696 variants loaded from bim file


# Plot of the heterozygosity rate vs FMISS (proportion of missing SNPs per individual)
# See R script: 2_2000HIV_GWAS_QC_plots_29nov2022.R

# Subset samples of Asian ancestry
# Options in effect 
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--keep 2000HIV_asians.txt
--make-bed 
--chr 1-22
--exclude high-LD-regions-hg19-GRCh37.txt 
--indep-pairwise 50 5 0.2 
--out asians_indepSNP
# 479555 of 665394 variants removed

# Extract the pruned SNPs
# Options in effect 
--bfile  all_ethnicities_mind0.025_maf0.01_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--extract asians_indepSNP.prune.in
--keep 2000HIV_asians.txt
--chr 1-22
--make-bed 
--out asians_only_autosomes_pruned
# 185839 variants and 85 people pass filters and QC

# Calculation of heterozygosity rate to identify individuals with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes
# Options in effect 
--bfile asians_only_autosomes_pruned
--het
--out asians_het_only_autosomes # Use this file to calculate heterozygosity rate
# Report written to asians_het_only_autosomes.het

#  samples of Asian ancestry were identified with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes

# Calculate missingness using only autosomes
# For all ethnicities
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--chr 1-22
--missing
--keep 2000HIV_asians.txt
--out asians_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
# 686171out of 714696 variants loaded from bim file


# Plot of the heterozygosity rate vs FMISS (proportion of missing SNPs per individual)
# See R script: 2_2000HIV_GWAS_QC_plots_29nov2022.R

# Subset samples of Hispanic ancestry
# Options in effect 
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--keep 2000HIV_hispanics.txt
--make-bed 
--chr 1-22
--exclude high-LD-regions-hg19-GRCh37.txt 
--indep-pairwise 50 5 0.2 
--out hispanics_indepSNP
# 513607 of 686171 variants removed

# Extract the pruned SNPs
# Options in effect 
--bfile  all_ethnicities_mind0.025_maf0.01_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--extract hispanics_indepSNP.prune.in
--keep 2000HIV_hispanics.txt
--chr 1-22
--make-bed 
--out hispanics_only_autosomes_pruned
# 172564 variants and 47 people pass filters and QC

# Calculation of heterozygosity rate to identify individuals with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes
# Options in effect 
--bfile hispanics_only_autosomes_pruned
--het
--out hispanics_het_only_autosomes # Use this file to calculate heterozygosity rate
# Report written to hispanics_het_only_autosomes.het

#  samples of Hispanic ancestry were identified with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes

# Calculate missingness using only autosomes
# For all ethnicities
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--chr 1-22
--missing
--keep 2000HIV_hispanics.txt
--out hispanics_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
# 686171 out of 714696 variants loaded from bim file


# Plot of the heterozygosity rate vs FMISS (proportion of missing SNPs per individual)
# See R script: 2_2000HIV_GWAS_QC_plots_29nov2022.R

# Subset samples of Native American ancestry
# Options in effect 
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--keep 2000HIV_native_american.txt
--make-bed 
--chr 1-22
--exclude high-LD-regions-hg19-GRCh37.txt 
--indep-pairwise 50 5 0.2 
--out native_americans_indepSNP
# 672932 of 686171 variants removed

# Extract the pruned SNPs
# Options in effect 
--bfile  all_ethnicities_mind0.025_maf0.01_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--extract native_americans_indepSNP.prune.in
--keep 2000HIV_native_american.txt
--chr 1-22
--make-bed 
--out native_americans_only_autosomes_pruned
# 13239 variants and 3 people pass filters and QC

# Calculation of heterozygosity rate to identify individuals with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes
# Options in effect 
--bfile native_americans_only_autosomes_pruned
--het
--out native_americans_het_only_autosomes # Use this file to calculate heterozygosity rate
# Report written to native_americans_het_only_autosomes.het

#  samples of Native American ancestry were identified with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes

# Calculate missingness using only autosomes
# For all ethnicities
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--chr 1-22
--missing
--keep 2000HIV_native_american.txt
--out native_americans_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
# 686171 out of 714696 variants loaded from bim file


# Plot of the heterozygosity rate vs FMISS (proportion of missing SNPs per individual)
# See R script: 2_2000HIV_GWAS_QC_plots_29nov2022.R

# Subset samples of Mixed ancestry
# Options in effect 
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--keep 2000HIV_mixed.txt
--make-bed 
--chr 1-22
--exclude high-LD-regions-hg19-GRCh37.txt 
--indep-pairwise 50 5 0.2 
--out mixed_indepSNP
# 408136 of 686171 variants removed

# Extract the pruned SNPs
# Options in effect 
--bfile  all_ethnicities_mind0.025_maf0.01_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 
--extract mixed_indepSNP.prune.in
--keep 2000HIV_mixed.txt
--chr 1-22
--make-bed 
--out mixed_only_autosomes_pruned
# 278035 variants and 128 people pass filters and QC

# Calculation of heterozygosity rate to identify individuals with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes
# Options in effect 
--bfile mixed_only_autosomes_pruned
--het
--out mixed_het_only_autosomes # Use this file to calculate heterozygosity rate
# Report written to mixed_het_only_autosomes.het

#  samples of Mixed ancestry were identified with excess heterozygosity (deviating more than 3 sd from the mean) using only the pruned SNPs from autosomes

# Calculate missingness using only autosomes
# For all ethnicities
--bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23
--chr 1-22
--missing
--keep 2000HIV_mixed.txt
--out mixed_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes
# 686171 out of 714696 variants loaded from bim file

# Plot of the heterozygosity rate vs FMISS (proportion of missing SNPs per individual)
# See R script: 2_2000HIV_GWAS_QC_plots_29nov2022.R

# Number of individuals with excess heterozygosity (deviating more than 3 sd from the mean)
White ancestry: 32 (whites_fail_het_qc_21dec2022.txt)
Black ancestry: 3 (blacks_fail_het_qc_21dec2022.txt)
Asian ancestry: 1 (asians_fail_het_qc_21dec2022.txt)
Hispanic ancestry: 1 (hispanics_fail_het_qc_21dec2022.txt)
Native american ancestry: 0
Mixed ancestry: 2 (mixed_fail_het_qc_21dec2022.txt)
# In total, 39 outliers with excess heterozygosity

#  Do PCA analysis - see script: 2000HIV_pca_analysis_all_ethnicities_with_1000G_hg37.sh

# Remove outliers based on heterozygosity rate (n = 39 samples) from the 2000HIV cohort of all ethncities
qcdir='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/pca_ancestry_1000G_all_ethnicities'

plink --bfile  $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 --output-chr chrM --remove all_het_outliers.txt --make-bed --out all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean

# 714696 variants and 1864 people pass filters and QC.

###############################################################################################################################################################################
# Lift the SNPs from GRCh37 to GRCh38 using hte liftOver tool, align againt the TopMed panel and submit for imputation using the TopMed server

# All ethnicities cohort ----------------------------------------------------------------------------------------------------------------------------------------------
qcdir='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data'

# Extract the SNP position and chromosome from the .bim file 
awk '{print$1,$4,$2}'  $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean.bim > all_ethnicities_clean_snps.bed
#  714696 SNPs in GRCh37 build

# Open R studio and import the bed file with the SNP chromorome and position
setwd("/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/")

library(data.table)
dir <- "/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/"
all_ethnicities_snps  <- fread(paste0(dir,"all_ethnicities_clean_snps.bed"), data.table = F) 
head(all_ethnicities_snps)
dim(all_ethnicities_snps)
# [1] 714696       3
tail(all_ethnicities_snps)

# Prepare .bed file as described here: https://genome.sph.umich.edu/wiki/LiftOver#Binary_liftOver_tool
all_ethnicities_snps$V4 <- all_ethnicities_snps$V2+1
head(all_ethnicities_snps)
all_ethnicities_final <- all_ethnicities_snps[,c("V1","V2","V4","V3")]
head(all_ethnicities_final)

dim(all_ethnicities_final)
# [1] 714696     4
tail(all_ethnicities_final)

# Convert your numbers to strings with formatting as you require, then using the argument quote = FALSE in the call to write.table.
all_ethnicities_final$V4 <- format(all_ethnicities_final$V4, scientific = FALSE)

# save
# write.table(all_ethnicities_final, file = "all_ethnicities_snps_final_GRCh37.bed", col.names = F, row.names = F, quote = F)

# Lift from hg19 to hg18
# hg19 to hg18
qcdir='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data'

./liftOver $qcdir/all_ethnicities_snps_final_GRCh37.bed  hg19ToHg18.over.chain.gz $qcdir/output_all_ethnicities_snps_final_hg18.bed $qcdir/unlifted_all_ethnicities_snps_final_hg18.bed


# Number of SNPs that were lifted from hg19 to hg18 in the output file (all_ethnicities cohort): 714399  SNPs
# Number of SNPs that were not lifted from hg19 to hg18: 297  SNPs

# Do second lift from hg18 to hg38
./liftOver $qcdir/output_all_ethnicities_snps_final_hg18.bed   hg18ToHg38.over.chain.gz  $qcdir/output2_all_ethnicities_snps_final_hg38.bed $qcdir/unlifted2_all_ethnicities_snps_final_hg38.bed

# Number of SNPs that were lifted from hg18 to hg38 in the output file (all_ethnicities cohort):   714292  SNPs (including 4 SNPs with invalid chromosome code; see below steps)
# Number of SNPs that were not lifted from hg18 to hg38:  107 SNPs
# Second column of the output file is the lifted SNPs to hg38 (GRCh38)

# Open R studio
# Convert lifted .bed file back to .map file
# Rearrange columns of lifted output.bed file to obtain.map file in the new build as follows:
# 1	GSA-rs116587930	0	727841
# 1	rs3131972	    0	752721

library(data.table)
library(stringr)

dir <- "/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/"
all_ethnicities_lifted_snps  <- fread(paste0(dir,"output2_all_ethnicities_snps_final_hg38.bed"), data.table = F) 
head(all_ethnicities_lifted_snps)
dim(all_ethnicities_lifted_snps)
# [1] 714292      4
tail(all_ethnicities_lifted_snps)

# Remove the "chr" from chr1-22
all_ethnicities_lifted_snps2 <- str_split_fixed(all_ethnicities_lifted_snps$V1, "chr", 2)
head(all_ethnicities_lifted_snps2)

all_ethnicities_lifted_snps3 <- cbind(all_ethnicities_lifted_snps2, all_ethnicities_lifted_snps)
head(all_ethnicities_lifted_snps3)

all_ethnicities_lifted_snps3$V5 <- '0'

# Most analyses do not require a genetic map to be specified in any case; 
# specifying a genetic (cM) map is most crucial for a set of analyses that look for shared segments between individuals. 
# For basic association testing, the genetic distance column can be set at 0.
# Above, we specify the genetic distance at 0.

# Remove columns 1 and 3
all_ethnicities_lifted_snps4 <- all_ethnicities_lifted_snps3[,c("2", "V4","V5","V2")]
dim(all_ethnicities_lifted_snps4)
# [1] 714292        4
tail(all_ethnicities_lifted_snps4)
head(all_ethnicities_lifted_snps4)

# Remove SNPs with invalide chromosome code
# chr4_GL000008v2_random	81626	81627	GSA-1:142617060
# chr4_GL000008v2_random	108295	108296	GSA-rs375872403
# chr14_GL000009v2_random	152859	152860	GSA-1:143495675
# chr14_GL000009v2_random	177760	177761	GSA-1:143520576

colnames(all_ethnicities_lifted_snps4)[1] <- 'chr'

all_ethnicities_lifted_snps5 <- all_ethnicities_lifted_snps4[!grepl("*_random", all_ethnicities_lifted_snps4$chr),]

dim(all_ethnicities_lifted_snps4 )
#  [1] 714292      4

dim(all_ethnicities_lifted_snps5)
# [1] 714288      4 # 4 SNPs were removed

# Save as txt file without the header
# write.table(all_ethnicities_lifted_snps5, file = "~/Desktop/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_toLift.map", col.names = F, row.names = F, quote = F, sep = "\t")

# As some genome positions cannot be lifted to the new version, we need to drop their corresponding columns from .ped file to keep consistency. 
# You can use PLINK --exclude those snps

# Extract unlifted SNPs using an awk-script
qcdir='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data'

awk '{print$4}' $qcdir/unlifted_all_ethnicities_snps_final_hg18.bed > $qcdir/all_ethnicities_unlifted_snps_hg18.txt
awk '{print$4}' $qcdir/unlifted2_all_ethnicities_snps_final_hg38.bed > $qcdir/all_ethnicities_unlifted_snps_hg38.txt


# Open R studio and import the file with the unlifted SNPs
dir <- "/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/"
all_ethnicities_unlifted_snps1  <- read.table(paste0(dir,"all_ethnicities_unlifted_snps_hg18.txt"), head = F) 
dim(all_ethnicities_unlifted_snps1)
# [1]  297     1

all_ethnicities_unlifted_snps2  <- read.table(paste0(dir,"all_ethnicities_unlifted_snps_hg38.txt"), head = F) 
dim(all_ethnicities_unlifted_snps2)
# [1]  107   1

# Add the unlifted SNPs in one file and save in a txt file
all_ethnicities_unlifted_snps_all <- rbind(all_ethnicities_unlifted_snps1, all_ethnicities_unlifted_snps2)
dim(all_ethnicities_unlifted_snps_all)
#[1] 404   1

# Add 4 more SNPs to the list with invalid chromosome code
# chr4_GL000008v2_random	81626	81627	GSA-1:142617060
# chr4_GL000008v2_random	108295	108296	GSA-rs375872403
# chr14_GL000009v2_random	152859	152860	GSA-1:143495675
# chr14_GL000009v2_random	177760	177761	GSA-1:143520576

# Note that there are some invalid chromosome codes ('chr*_random') in GRCh38 build. In total, it is 4 SNPs. 
# We will remove them when preparing the .map file in R studio (see above)


invalid_chr <- data.frame(V1=c("GSA-1:142617060","GSA-rs375872403","GSA-1:143495675","GSA-1:143520576"))
all_ethnicities_unlifted_snps_all_final <- rbind(all_ethnicities_unlifted_snps_all, invalid_chr)
dim(all_ethnicities_unlifted_snps_all_final)
# [1] 408   1

# write.table(all_ethnicities_unlifted_snps_all_final, file = "all_ethnicities_unlifted_snps_all_to_exclude.txt", row.names = F, col.names = F, sep = "\t", quote = F)

# Remove the unlifted SNPs and convert to ped files
plink --bfile $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean --exclude $qcdir/all_ethnicities_unlifted_snps_all_to_exclude.txt --recode --tab --out $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_toLift
# 714696 variants loaded from .bim file.
# 714288 variants and 1864 people pass filters and QC.

# Use the lifted .map file and convert to bed file
plink --file $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_toLift  --make-bed --out $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38 
# 714288 variants and 1864 people pass filters and QC.


# Similar to the human reference build, dbSNP also have different versions. You may consider change rs number from the old dbSNP version to 
# new dbSNP version depending on your needs. 

# Files to use for alingment to TopMed
all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38.bed
all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38.bim
all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38.fam

####################################################################################################################################################################################
####################################################################################################################################################################################

## Below, we aligned the pre-imputed genetic data to the TopMed reference panel and then submit for imputation using the TopMed reference panel -------------------------------------------------------------------------------------------------------- 
# Step 14 - Calculate allele frequencies before checking the QC'd file against HRC reference genome in advance of imputation
# all_ethnicities cohort
plink --bfile all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38  --freq --out all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38.freq

#  Align the pre-imputed genetic data to the TopMed panel (chromosomes should be coded as: 1-23)
#  Script to check plink .bim files against TopMed reference genome for strand, id names, positions, alleles, ref/alt assignment (source: https://www.well.ox.ac.uk/~wrayner/tools/)

# Once the TOPMed reference panel downloaded, the VCF was converted to an HRC formatted reference legend using the code here: CreateTOPMed.zip (see instructions here: https://www.well.ox.ac.uk/~wrayner/tools/)
# Usage: ./CreateTOPMed.pl -i ALL.TOPMed_freeze5_hg38_dbSNP.vcf.gz
# By default this will create a file filtered for variants flagged as PASS only, if you wish to use all variants the -a flag overrides this.
# To override the default output file naming use -o filename.

# Files required in the same directory to align:
# HRC-1000G-check-bim.pl
# PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab

qcdir='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data'

perl HRC-1000G-check-bim.pl -b $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38.bim -f $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38.freq.frq -h -r PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab -o TopMed_all_ethnicities_filtered_2000HIV_GSA796_chr1to23.clean_GRCh38

# Freq threshold: default: 0.2; Frequency difference to use when checking allele frequency of data set versus reference 
# Note, that I used the defaul output file naming (so we read HRC), and did not use -o filename to override the default output file naming.


# Note if your input data is GRCh38hg38  ensure chromosomes are encoded with prefix 'chr' (e.g. chr20).      https://imputationserver.readthedocs.io
# Edit the 'Run-plink.sh' and add flag '--output-chr chrM' (from line 5 onwords) so chromosomes are encoded with prefix 'chr'; use script 'Run-plink_for_TopMed_server.sh'

# Move the plink and vcf files to separate folders
qcdir2='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/all_ethnicities_HRC-1000G-check-bim-v4.3.0_jan2023/TopMed_all_ethnicities_filtered_2000HIV_GSA796_chr1to23.clean_GRCh38'

mkdir $qcdir2/plink
mkdir $qcdir2/vcf

# plink files
for chr in {1..23}
do
mv $qcdir2/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr${chr}.bim  $qcdir2/plink/
done

for chr in {1..23}
do
mv $qcdir2/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr${chr}.fam  $qcdir2/plink/
done

for chr in {1..23}
do
mv $qcdir2/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr${chr}.log  $qcdir2/plink/
done

for chr in {1..23}
do
mv $qcdir2/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr${chr}.bed  $qcdir2/plink/
done

for chr in {1..23}
do
mv $qcdir2/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr${chr}.nosex  $qcdir2/plink/
done


mv $qcdir2/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr23.hh  $qcdir2/plink/


# vcf files
for chr in {1..23}
do
mv $qcdir2/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr${chr}.vcf  $qcdir2/vcf/
done

#  Create a sorted vcf.vcf.gz files 
# If using BCF tools, see (source: https://topmedimpute.readthedocs.io/en/latest/prepare-your-data.html)

inputFolder="/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/all_ethnicities_HRC-1000G-check-bim-v4.3.0_jan2023/TopMed_all_ethnicities_filtered_2000HIV_GSA796_chr1to23.clean_GRCh38/vcf/"
mkdir ${inputFolder}for_TopMed_imputation
outputFolder="${inputFolder}for_TopMed_imputation/"

# For all chromosomes
for chr in {1..23}
do
vcf-sort ${inputFolder}all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.clean_GRCh38-updated-chr${chr}.vcf | bgzip -c > ${outputFolder}chr_${chr}_GRCh38_2000HIV_all_ethnicities_clean_final.vcf.gz
done

# Submit pre-imputed genetic data for imputation using the Top imputation server (Minimac4) 1.6.6 - https://imputation.biodatacatalyst.nhlbi.nih.gov/#!

# Check for data preparation: https://topmedimpute.readthedocs.io/en/latest/prepare-your-data.html

# Imputation settings
# Job name: 3Jan2023 - 2000HIV cohort of all ancestries
# Reference panel: TOPMed r2
# Array Build: GRCh38/hg38
# Phasing: Eagle v2.4 (phased output)
# Population: vs. TOPMed Panel
# Mode: Quality Control & Imputation

# Input validation (Path: )
# 23 valid VCF file(s) found.
# chr_{1-23}_GRCh38_2000HIV_all_ethnicities_clean_final.vcf.gz

# Samples: 1864
# Chromosomes: 10 11 12 13 14 15 16 17 18 19 1 20 21 22 X 2 3 4 5 6 7 8 9
# SNPs: 582404
# Chunks: 309
# Datatype: unphased
# Build: hg38
# Reference Panel: apps@topmed-r2@1.0.0 (hg38)
# Population: all
# Phasing: eagle
# Mode: imputation

# Input Validation
# 23 valid VCF file(s) found.

# Samples: 1864
# Chromosomes: 10 11 12 13 14 15 16 17 18 19 1 20 21 22 X 2 3 4 5 6 7 8 9
# SNPs: 582404
# Chunks: 309
# Datatype: unphased
# Build: hg38
# Reference Panel: apps@topmed-r2@1.0.0 (hg38)
# Population: all
# Phasing: eagle
# Mode: imputation

# Quality Control
# Calculating QC Statistics

# Statistics:
# Alternative allele frequency > 0.5 sites: 81,862
# Reference Overlap: 90.87 %
# Match: 529,160
# Allele switch: 0
# Strand flip: 0
# Strand flip and allele switch: 0
# A/T, C/G genotypes: 0
# Filtered sites:
# Filter flag set: 0
# Invalid alleles: 0
# Multiallelic sites: 0
# Duplicated sites: 0
# NonSNP sites: 0
# Monomorphic sites: 0
# Allele mismatch: 57
# SNPs call rate < 90%: 0

# Excluded sites in total: 57
# Remaining sites in total: 529,160
# See snps-excluded.txt for details
# Typed only sites: 53,187
# See typed-only.txt for details

# Warning: 1 Chunk(s) excluded: < 3 SNPs (see chunks-excluded.txt for details).
# Remaining chunk(s): 308

# THE END ############################################################################################################################################################################


























