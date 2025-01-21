# Readme
# Created by: Vasiliki Matzaraki
# Created on: 2 January 2023

# Source code: https://meyer-lab-cshl.github.io/plinkQC/articles/AncestryCheck.html


# Ancestry estimation of the 2000HIV cohort, including all ethnicities, using 1000Genomes reference genome from build 37
# The identification of individuals of divergent ancestry can be achieved by combining the genotypes of the study population with genotypes of a reference dataset consisting of individuals from known ethnicities (1000Genomes, build 37)


# Set-up ------------------------------------------------------------------------------------------------------------------------

qcdir='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/pca_ancestry_1000G_all_ethnicities'
qcdir2='/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data'
refdir='/Users/vickymatzaraki/Desktop/Bioinformatics/public_data/genome_references/1000Genomes_phase3_hg37_2504samples'
# mkdir $qcdir/pca_ancestry_1000G_all_ethnicities

# Match study genotypes and reference data
# In order to compute joint principal components of the reference and study population, we’ll need to combine the two datasets. 
# The plink –merge function enables this merge, but requires the variants in the datasets to be matching by chromosome, position and alleles. 
# The following sections show how to extract the relevant data from the reference and study dataset and how to filter matching variants.

# Filter reference and study data for non A-T or G-C SNPs ------------------------------------------------------------------------------------------------------------------------
# We will use an awk script to find A→T and C→G SNPs. As these SNPs are more difficult to align and only a subset of SNPs is required for the analysis, we will remove them from 
# both the reference and study data set.

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23.bim  > \
    $qcdir/all_ethnicities_ac_gt_snps.txt

    # 6583 SNPs with GC/CG and AT/TA alleles

   
plink --bfile  $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23 \
      --exclude $qcdir/all_ethnicities_ac_gt_snps.txt \
      --make-bed \
      --out $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps

# Prune study data ---------------------------------------------------------------------------------------------------------------------------------------------------
# We will conduct principle component analysis on genetic variants that are pruned for variants in linkage disequilibrium (LD) with an 𝑟2>0.2 in a 50kb window. 
# The LD-pruned dataset is generated below, 
# using plink –indep-pairwise to compute the LD-variants; additionally exclude range is used to remove genomic ranges of known high-LD structure.
# The file with the genomic ranges of known high-LD structure was provided by the following link:
# https://github.com/meyer-lab-cshl/plinkQC/tree/master/inst/extdata/hihg-LD-regions-hg19-GRCh37.txt

plink --bfile  $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps \
      --exclude range  /Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/pca_ancestry_1000G_discovery/high-LD-regions-hg19-GRCh37.txt  \
      --indep-pairwise 50 2 0.2 \
      --out $qcdir/pruning1_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps



plink --bfile  $qcdir/all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps \
      --extract $qcdir/pruning1_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/LDpruned_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps



# Filter reference data for the same SNP set as in study
# We will use the list of pruned variants from the study sample to reduce the reference dataset to the size of the study samples:

plink --bfile  $qcdir2/pca_ancestry_1000G_discovery/maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps \
      --extract $qcdir/pruning1_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps.prune.in \
      --make-bed \
      --allow-extra-chr \
      --out $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps


# Check and correct chromosome mismatch
# The following section uses an awk-script to check that the variant IDs of the reference data have the same chromosome ID as the study data. 
# For computing the genetic PC, the annotation is not important, however, merging the files via PLINK will only work for variants with perfectly matching attributes.
# For simplicity, we update the pruned reference dataset. Note, that sex chromosomes are often encoded differently and might make the matching more difficult. 
# Again, for simplicity and since not crucial to the final task, we will ignore XY-encoded sex chromosomes (via sed -n '/^[XY]/!p').

 #awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
 #   ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
 #    $qcdir/pca_ancestry_1000G_all_ethnicities/LDpruned_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps.bim $qcdir/pca_ancestry_1000G_all_ethnicities/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.bim | \
 #    sed -n '/^[23]/!p' > $qcdir/pca_ancestry_1000G_all_ethnicities/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.pruned_toUpdateChr_v2.txt


# Include X chromosome
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/LDpruned_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps.bim $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.bim > $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.pruned_toUpdateChr.txt


# The variant IDs of the reference data (1000Genomes) have the same chromosome ID as the study data (2000HIV all_ethnicities cohort)

#plink --bfile $qcdir/pca_ancestry_1000G/$refname.pruned \
#      --update-chr $qcdir/$refname.toUpdateChr 1 2 \
#      --make-bed \
#      --out $qcdir/$refname.updateChr
# mv $qcdir/$refname.updateChr.log $qcdir/plink_log/$refname.updateChr.log


# Position mismatch
# Similar to the chromosome matching, we use an awk-script to find variants with mis-matching chromosomal positions.

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/LDpruned_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps.bim  $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.bim > \
    $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps_toUpdatePos.txt


# Possible allele flips
# Unlike chromosomal and base-pair annotation, mismatching allele-annotations will not only prevent the plink –merge, 
# but also mean that it is likely that actually a different genotype was measured. 
# Initially, we can use the following awk-script to check if non-matching allele codes are a simple case of allele flips.

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/LDpruned_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps.bim  $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.bim > \
    $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps_toFlip.txt


# Upate positions and flip alleles
# We use plink to update the mismatching positions and possible allele-flips identified above.

plink --bfile $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps \
      --update-map $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps_toUpdatePos.txt 1 2 \
      --flip $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps_toFlip.txt \
      --make-bed \
      --out $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.flipped

# 12191 SNPs flipped.

# Remove mismatches
# Any alleles that do not match after allele flipping, are identified and removed from the reference dataset.

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $qcdir/LDpruned_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps.bim $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.flipped.bim > \
    $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.flipped.mismatch.txt

plink --bfile $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.flipped\
      --exclude $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.flipped.mismatch.txt \
      --make-bed \
      --out $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.pruned.clean

# Merge study genotypes and reference data
# The matching study and reference dataset can now be merged into a combined dataset with plink –bmerge. 
# If all steps outlined above were conducted successfully, no mismatch errors should occur.
plink --bfile $qcdir/LDpruned_all_ethnicities_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_chr1to23_no_ac_gt_snps  \
      --bmerge $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.pruned.clean.bed $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.pruned.clean.bim \
         $qcdir/LDpruned_maf0.01_no_relatives_all_phase3_hg37_no_ac_gt_snps.pruned.clean.fam  \
      --make-bed \
      --out $qcdir/2000HIV_all_ethnicities.merge.1000G_build_GRCh37

# PCA on the merged data
# We can now run principal component analysis on the combined dataset using plink –pca which returns a .eigenvec file 
# with the family and individual ID in columns 1 and 2, followed by the first 20 principal components.


plink --bfile $qcdir/2000HIV_all_ethnicities.merge.1000G_build_GRCh37 \
      --pca 6 \
      --out $qcdir/pcas6_2000HIV_all_ethnicities.merge.1000G_build_GRCh37


# Prapare MDS plot using R studio - The following steps were performed in R studio -------------------------------------------------------------------------------------
# Import libraries
library(data.table)
library(ggplot2)
library(gridExtra)
library(lemon)
library(plyr)
library(tidyr)

# Set working directory
dir <- "/Users/vickymatzaraki/2000HIV_project/2000HIV_data_analysis/2000HIV_QC_and_imputation/final_QC_of_preimputed_genetic_data/pca_ancestry_1000G_all_ethnicities/"
dir_pop <- "/Users/vickymatzaraki/Desktop/Bioinformatics/public_data/genome_references/1000Genomes_phase3_hg37_2504samples/"
dir_pheno <- "/Users/vickymatzaraki/Data_cohorts_fromJan2021/2000HIV_study/2000HIV/ECFR_GUIDE_Baseline_parameters_version1.0/"

# Import the input with the principal components
pcaFile <- fread(paste0(dir,"pcas6_2000HIV_all_ethnicities.merge.1000G_build_GRCh37.eigenvec"), data.table = FALSE)
dim(pcaFile)
# [1] 4393    8
head(pcaFile)

colnames(pcaFile) <- c("FID", "IID", "PC1", "PC2",  "PC3",  "PC4",  "PC5", "PC6")
colnames(pcaFile) 

# Import the input with the 1000Genomes subpopulations
subpopulations <- fread(paste0(dir_pop,"igsr_1000_genomes_phase3_release.txt"), data.table = FALSE)
colnames(subpopulations)
# [1] "Sample name"           "Sex"                   "Biosample ID"          "Population code"       "Population name"      
# [6] "Superpopulation code"  "Superpopulation name"  "Population elastic ID" "Data collections"     
colnames(subpopulations)[1] <- "IID"

subpopulations_only <- subpopulations[,c("IID","Superpopulation code"),]

# Join
pcaFile_with_pops <- join(pcaFile, subpopulations_only, by = "IID")
dim(pcaFile_with_pops )
# [1] 4393   9

# Import the phenotype information
ethnicity <- fread(paste0(dir_pheno,"221110_2000hiv_study_export_processed_2.0_SIMPLIFIED.csv"), data.table = FALSE)

# [1] "ID"                                 "CENTER"                             "COHORT"                            
# [4] "SEX_BIRTH"                          "GEN"                                "AGE"                               
# [7] "DATE_VISIT"                         "PANDEMIC_BEFOREAFTER"               "ETHNICITY"

ethnicity <- ethnicity[,c("ID","ETHNICITY")]
# Change column names
colnames(ethnicity)[1] <- 'IID'

# Extract ethnicity for the 2000HIV cohort
pcaFile_with_pops_v2 <- join(pcaFile_with_pops, ethnicity, by = "IID")
head(pcaFile_with_pops_v2)

# Combine columns (`Superpopulation code` and "ETHNICITY") to remove NAs so ethnicity of the 2000HVI and 1000G is in one column (new)
pcaFile_with_pops_v2$Superpopulation =  pcaFile_with_pops_v2$`Superpopulation code`
pcaFile_with_pops_v2$Superpopulation[!is.na(pcaFile_with_pops_v2$ETHNICITY)] = pcaFile_with_pops_v2$ETHNICITY[!is.na(pcaFile_with_pops_v2$ETHNICITY)]  # merge with y

unique(pcaFile_with_pops_v2$Superpopulation)
# [1] "EUR"             "EAS"             "AMR"             "SAS"             "EUR,AFR"         "AFR"            
# [7] "White"           "Black"           ""                "Asian"           NA                "Mixed"          
# [13] "Hispanic"        "Native American"

# Count the samples per ethnicity group
white <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "White"))
length(white$Superpopulation)
# [1] 1374

black <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Black"))
length(black$Superpopulation)
#[1] 184

asian <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Asian"))
length(asian$Superpopulation)
#[1] 85

mixed <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Mixed"))
length(mixed$Superpopulation)
# [1] 128

hispanic <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Hispanic"))
length(hispanic$Superpopulation)
# [1] 47
nativeamerican <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Native American"))
length(nativeamerican$Superpopulation)
# [1] 3

# Prepare PCA plots per ethnicity
unique(pcaFile_with_pops_v2$Superpopulation)
# White ancestry
white_plot <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "White" | pcaFile_with_pops_v2$Superpopulation == "EUR" |  pcaFile_with_pops_v2$Superpopulation ==  "EAS" |  pcaFile_with_pops_v2$Superpopulation == "AMR" | pcaFile_with_pops_v2$Superpopulation == "SAS" | pcaFile_with_pops_v2$Superpopulation == "EUR,AFR" | pcaFile_with_pops_v2$Superpopulation == "AFR"))

legend_title <- "Population"
p1 = ggplot(white_plot, aes(PC1, PC2,  shape = Superpopulation, color = Superpopulation, group = Superpopulation)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()+
  ggtitle("2000HIV cohort of European ancestry vs 1000 Genomes (GRCh37)")

p2 = ggplot(white_plot, aes(PC2, PC3,  shape = Superpopulation, color = Superpopulation, group = Superpopulation)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

p3 = ggplot(white_plot, aes(PC3, PC4,  shape = Superpopulation, color = Superpopulation, group = Superpopulation)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

# Combine the plots
nt <- theme(legend.position='none')
p4 <- grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)

# pdf(file = "~/Desktop/PCA_plot_2000HIV_all_White_ancestry_vs_1000Genomes_GRCh37_2jan2022.pdf", width = 8, height =8 )
grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)
# dev.off()


# black ancestry
black_plot <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Black" | pcaFile_with_pops_v2$Superpopulation == "EUR" |  pcaFile_with_pops_v2$Superpopulation ==  "EAS" |  pcaFile_with_pops_v2$Superpopulation == "AMR" | pcaFile_with_pops_v2$Superpopulation == "SAS" | pcaFile_with_pops_v2$Superpopulation == "EUR,AFR" | pcaFile_with_pops_v2$Superpopulation == "AFR"))

legend_title <- "Population"

black_plot$Population <- factor(black_plot$Superpopulation, levels = c("AFR", "AMR", "EAS", "EUR", "EUR,AFR", "SAS","Black"))

p1 = ggplot(black_plot, aes(PC1, PC2,  shape = black_plot$Population, color = black_plot$Population, group = black_plot$Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()+
  ggtitle("2000HIV cohort of African ancestry vs 1000 Genomes (GRCh37)")+
  scale_fill_discrete(name = "Population")

p2 = ggplot(black_plot, aes(PC2, PC3,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

p3 = ggplot(black_plot, aes(PC3, PC4,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

# Combine the plots
nt <- theme(legend.position='none')
p4 <- grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)

#pdf(file = "~/Desktop/PCA_plot_2000HIV_all_Black_ancestry_vs_1000Genomes_GRCh37_2jan2022.pdf", width = 8, height =8 )
grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)
#dev.off()

# asian ancestry
asian_plot <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Asian" | pcaFile_with_pops_v2$Superpopulation == "EUR" |  pcaFile_with_pops_v2$Superpopulation ==  "EAS" |  pcaFile_with_pops_v2$Superpopulation == "AMR" | pcaFile_with_pops_v2$Superpopulation == "SAS" | pcaFile_with_pops_v2$Superpopulation == "EUR,AFR" | pcaFile_with_pops_v2$Superpopulation == "AFR"))

legend_title <- "Population"

asian_plot$Population <- factor(asian_plot$Superpopulation, levels = c("AFR", "AMR", "EAS", "EUR", "EUR,AFR", "SAS","Asian"))

p1 = ggplot(asian_plot, aes(PC1, PC2,  shape = asian_plot$Population, color = asian_plot$Population, group = asian_plot$Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()+
  ggtitle("2000HIV cohort of Asian ancestry vs 1000 Genomes (GRCh37)")+
  scale_fill_discrete(name = "Population")

p2 = ggplot(asian_plot, aes(PC2, PC3,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

p3 = ggplot(asian_plot, aes(PC3, PC4,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

# Combine the plots
nt <- theme(legend.position='none')
p4 <- grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)

#pdf(file = "~/Desktop/PCA_plot_2000HIV_all_Asian_ancestry_vs_1000Genomes_GRCh37_2jan2022.pdf", width = 8, height =8 )
grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)
#dev.off()


# mixed ancestry
mixed_plot <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Mixed" | pcaFile_with_pops_v2$Superpopulation == "EUR" |  pcaFile_with_pops_v2$Superpopulation ==  "EAS" |  pcaFile_with_pops_v2$Superpopulation == "AMR" | pcaFile_with_pops_v2$Superpopulation == "SAS" | pcaFile_with_pops_v2$Superpopulation == "EUR,AFR" | pcaFile_with_pops_v2$Superpopulation == "AFR"))

legend_title <- "Population"

mixed_plot$Population <- factor(mixed_plot$Superpopulation, levels = c("AFR", "AMR", "EAS", "EUR", "EUR,AFR", "SAS","Mixed"))

p1 = ggplot(mixed_plot, aes(PC1, PC2,  shape = mixed_plot$Population, color = mixed_plot$Population, group = mixed_plot$Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()+
  ggtitle("2000HIV cohort of Mixed ancestry vs 1000 Genomes (GRCh37)")+
  scale_fill_discrete(name = "Population")

p2 = ggplot(mixed_plot, aes(PC2, PC3,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

p3 = ggplot(mixed_plot, aes(PC3, PC4,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

# Combine the plots
nt <- theme(legend.position='none')
p4 <- grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)

#pdf(file = "~/Desktop/PCA_plot_2000HIV_all_Mixed_ancestry_vs_1000Genomes_GRCh37_2jan2022.pdf", width = 8, height =8 )
grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)
#dev.off()


# hispanic ancestry
hispanic_plot <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Hispanic" | pcaFile_with_pops_v2$Superpopulation == "EUR" |  pcaFile_with_pops_v2$Superpopulation ==  "EAS" |  pcaFile_with_pops_v2$Superpopulation == "AMR" | pcaFile_with_pops_v2$Superpopulation == "SAS" | pcaFile_with_pops_v2$Superpopulation == "EUR,AFR" | pcaFile_with_pops_v2$Superpopulation == "AFR"))

legend_title <- "Population"

hispanic_plot$Population <- factor(hispanic_plot$Superpopulation, levels = c("AFR", "AMR", "EAS", "EUR", "EUR,AFR", "SAS","Hispanic"))

p1 = ggplot(hispanic_plot, aes(PC1, PC2,  shape = hispanic_plot$Population, color = hispanic_plot$Population, group = hispanic_plot$Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()+
  ggtitle("2000HIV cohort of Hispanic ancestry vs 1000 Genomes (GRCh37)")+
  scale_fill_discrete(name = "Population")

p2 = ggplot(hispanic_plot, aes(PC2, PC3,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

p3 = ggplot(hispanic_plot, aes(PC3, PC4,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

# Combine the plots
nt <- theme(legend.position='none')
p4 <- grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)

#pdf(file = "~/Desktop/PCA_plot_2000HIV_all_Hispanic_ancestry_vs_1000Genomes_GRCh37_2jan2022.pdf", width = 8, height =8 )
grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)
#dev.off()


# Native american ancestry
nativeamerican_plot <- subset(pcaFile_with_pops_v2, (pcaFile_with_pops_v2$Superpopulation == "Native American" | pcaFile_with_pops_v2$Superpopulation == "EUR" |  pcaFile_with_pops_v2$Superpopulation ==  "EAS" |  pcaFile_with_pops_v2$Superpopulation == "AMR" | pcaFile_with_pops_v2$Superpopulation == "SAS" | pcaFile_with_pops_v2$Superpopulation == "EUR,AFR" | pcaFile_with_pops_v2$Superpopulation == "AFR"))

legend_title <- "Population"

nativeamerican_plot$Population <- factor(nativeamerican_plot$Superpopulation, levels = c("AFR", "AMR", "EAS", "EUR", "EUR,AFR", "SAS","Native American"))

p1 = ggplot(nativeamerican_plot, aes(PC1, PC2,  shape = nativeamerican_plot$Population, color = nativeamerican_plot$Population, group = nativeamerican_plot$Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()+
  ggtitle("2000HIV cohort of Native American ancestry vs 1000 Genomes (GRCh37)")+
  scale_fill_discrete(name = "Population")

p2 = ggplot(nativeamerican_plot, aes(PC2, PC3,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

p3 = ggplot(nativeamerican_plot, aes(PC3, PC4,  shape = Population, color = Population, group = Population)) +
  geom_point()+
  scale_colour_manual(
    values = c("red", "yellow", "black", "orange","blue","purple","darkgreen")) +   
  scale_shape_manual( values = c(19, 19,19,19,19,19,17))+
  theme_bw()

# Combine the plots
nt <- theme(legend.position='none')
p4 <- grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)

#pdf(file = "~/Desktop/PCA_plot_2000HIV_all_Native_American_ancestry_vs_1000Genomes_GRCh37_2jan2022.pdf", width = 8, height =8 )
grid_arrange_shared_legend(p1,p2,p3, ncol=1, nrow=3)
#dev.off()

# THE END ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


