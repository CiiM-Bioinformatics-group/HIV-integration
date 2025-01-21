# Readme
# Created on: 29 November 2022
# Source of raw genotype data: /mnt/workspace/1. Original data/1.Raw/2000HIV/2000HIV_GWAS_GSAData/GSA796/PLINK_220222_0252
# File names: GSA796.phenotype, GSA796.ped, GSA796.map, GSA796.bat, GSA796.script

# CREATE PLOTS DURING QC OF THE PRE-IMPUTED GENETIC DATA AND OTHER DATA HANDLING DURING QC ------------------------------------------------------------------------------------------------------------------------------------------------
# Load libraries ------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(plyr)
library(ggrepel)
library(ggplot2)
library(dplyr)

# Set working directory
setwd("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/")
dir <- "Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/"
dir_plots <- "Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots"

# Generate plots to visualize the missingness results per cohort-------------------------------------------------------------------------------------
# Upload the .imiss and lmiss files into R
indmiss_all <- fread(paste0(dir,"miss_2000HIV_GSA796_all_chromosomes.imiss"), data.table = F) 
head(indmiss_all)

snpmiss_all <- fread(paste0(dir,"miss_2000HIV_GSA796_all_chromosomes.lmiss"), data.table = F)
head(snpmiss_all)

#pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/histimiss_individual_missingenss_for_all_ethnicities.pdf") #indicates pdf format and gives title to file
hist(indmiss_all[,6],main="Histogram individual missingness, 2000HIV cohort (n = 1914)") #selects column 6, names header of file
#dev.off()

#pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/histimiss_SNP_missingenss_for_all_ethnicities.pdf") 
hist(snpmiss_all[,5],main="Histogram SNP missingness, 2000HIV (n = 1914)")  
#dev.off() # shuts down the current device


# Extract sample ids for only females from the all_ethnicities cohort
# we coded females with 2  and males with 1 
fam_all_ethnicities <- fread(paste0(dir, "all_ethnicities_geno0.05_european_updated_sex_new2_2000HIV_GSA796_only_autosomes.fam"), data.table = F)
dim(fam_all_ethnicities)
# 1914    6
head(fam_all_ethnicities)

fam_all_ethnicities_females <- subset(fam_all_ethnicities, fam_all_ethnicities$V5 == "2")
dim(fam_all_ethnicities_females)
# [1] 307   6

fam_all_ethnicities_females_final <- fam_all_ethnicities_females[,c("V1","V2"),]
dim(fam_all_ethnicities_females_final)
head(fam_all_ethnicities_females_final)

# save 
# write.table(fam_all_ethnicities_females_final, file = "only_females_for_all_ethnicities_cohort_1dec2022.txt", quote = F, row.names = F, col.names = F)

# Extract SNPs that removed due to Hardy-Weinberg exact test per cohort --------------------------------------------------------------------------------------------------------
# All ethnicities - autosomes
# Extract SNPs that removed due to Hardy-Weinberg exact test
snps_all_ethnicities_before_autosomes <- fread(paste0(dir, "all_ethnicities_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes.bim"), data.table = F)
head(snps_all_ethnicities_before_autosomes)
dim(snps_all_ethnicities_before_autosomes)
# [1] 686708      3

# Identify autosomal SNPs that failed HWE in samples of white ancestry
snps_whites_after_autosomes <- fread(paste0(dir, "whites_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_gsa796_only_autosomes.bim"),data.table = F)
head(snps_whites_after_autosomes)
dim(snps_whites_after_autosomes)
# [1] 686237      6

snps_failed_hwe_whites_autosomes <- anti_join(snps_all_ethnicities_before_autosomes, snps_whites_after_autosomes, by = "V2")
dim(snps_failed_hwe_whites_autosomes)
# [1] 471    6 # 471 variants failed HWE*

snps_failed_hwe_whites_autosomes_final <- as.data.frame(snps_failed_hwe_whites_autosomes[,c(2)])
head(snps_failed_hwe_whites_autosomes_final)
dim(snps_failed_hwe_whites_autosomes_final)
# [1] 471     1

colnames(snps_failed_hwe_whites_autosomes_final) <- 'SNP'

# Identify autosomal SNPs that failed HWE in samples of black ancestry
snps_blacks_after_autosomes <- fread(paste0(dir, "blacks_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_gsa796_only_autosomes.bim"),data.table = F)
head(snps_blacks_after_autosomes)
dim(snps_blacks_after_autosomes)
# [1] 686658      6

snps_failed_hwe_blacks_autosomes <- anti_join(snps_all_ethnicities_before_autosomes, snps_blacks_after_autosomes, by = "V2")
dim(snps_failed_hwe_blacks_autosomes)
# [1] 50    6 # 50 variants failed HWE*

snps_failed_hwe_blacks_autosomes_final <- as.data.frame(snps_failed_hwe_blacks_autosomes[,c(2)])
head(snps_failed_hwe_blacks_autosomes_final)
dim(snps_failed_hwe_blacks_autosomes_final)
# [1] 50     1

colnames(snps_failed_hwe_blacks_autosomes_final) <- 'SNP'

# Identify autosomal SNPs that failed HWE in samples of Asian ancestry
snps_asians_after_autosomes <- fread(paste0(dir, "asians_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_gsa796_only_autosomes.bim"),data.table = F)
head(snps_asians_after_autosomes)
dim(snps_asians_after_autosomes)
# [1] 686697       6

snps_failed_hwe_asians_autosomes <- anti_join(snps_all_ethnicities_before_autosomes, snps_asians_after_autosomes, by = "V2")
dim(snps_failed_hwe_asians_autosomes)
# [1] 11    6 # 11 variants failed HWE*

snps_failed_hwe_asians_autosomes_final <- as.data.frame(snps_failed_hwe_asians_autosomes[,c(2)])
head(snps_failed_hwe_asians_autosomes_final)
dim(snps_failed_hwe_asians_autosomes_final)
# [1] 11     1

colnames(snps_failed_hwe_asians_autosomes_final) <- 'SNP'

# Identify autosomal SNPs that failed HWE in samples of Hispanic ancestry
snps_hispanics_after_autosomes <- fread(paste0(dir, "hispanics_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_gsa796_only_autosomes.bim"),data.table = F)
head(snps_hispanics_after_autosomes)
dim(snps_hispanics_after_autosomes)
# [1] 686706       6

snps_failed_hwe_hispanics_autosomes <- anti_join(snps_all_ethnicities_before_autosomes, snps_hispanics_after_autosomes, by = "V2")
dim(snps_failed_hwe_hispanics_autosomes)
# [1] 2    6 # 2 variants failed HWE*

snps_failed_hwe_hispanics_autosomes_final <- as.data.frame(snps_failed_hwe_hispanics_autosomes[,c(2)])
head(snps_failed_hwe_hispanics_autosomes_final)
dim(snps_failed_hwe_hispanics_autosomes_final)
# [1] 2     1

colnames(snps_failed_hwe_hispanics_autosomes_final) <- 'SNP'

# Identify autosomal SNPs that failed HWE in samples of Mixed ancestry
snps_mixed_after_autosomes <- fread(paste0(dir, "mixed_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_gsa796_only_autosomes.bim"),data.table = F)
head(snps_mixed_after_autosomes)
dim(snps_mixed_after_autosomes)
# [1] 686689        6

snps_failed_hwe_mixed_autosomes <- anti_join(snps_all_ethnicities_before_autosomes, snps_mixed_after_autosomes, by = "V2")
dim(snps_failed_hwe_mixed_autosomes)
# [1] 19    6 # 19 variants failed HWE*

snps_failed_hwe_mixed_autosomes_final <- as.data.frame(snps_failed_hwe_mixed_autosomes[,c(2)])
head(snps_failed_hwe_mixed_autosomes_final)
dim(snps_failed_hwe_mixed_autosomes_final)
# [1] 19     1

colnames(snps_failed_hwe_mixed_autosomes_final) <- 'SNP'

# All ethnicities - X_chromosome
# Extract SNPs that removed due to Hardy-Weinberg exact test
snps_all_ethnicities_before_X_chrom <- fread(paste0(dir, "all_ethnicities_females_geno0.05_updated_sex_new2_2000HIV_GSA796_X_chromosome.bim"), data.table = F)
head(snps_all_ethnicities_before_X_chrom)
dim(snps_all_ethnicities_before_X_chrom)
# [1] 28526     3

snps_blacks_after_X_chrom <- fread(paste0(dir, "blacks_females_hwe1minus6_geno0.05_european_updated_sex_new2_2000HIV_GSA796_X_chromosome.bim"),data.table = F)
head(snps_blacks_after_X_chrom)
dim(snps_blacks_after_X_chrom)
# [1] 28525      6

snps_failed_hwe_blacks_X_chrom <- anti_join(snps_all_ethnicities_before_X_chrom, snps_blacks_after_X_chrom, by = "V2")
dim(snps_failed_hwe_blacks_X_chrom)
# [1]  1     6 #  1  variants failed HWE*

snps_failed_hwe_blacks_X_chrom_final <- as.data.frame(snps_failed_hwe_blacks_X_chrom[,c(2)])
head(snps_failed_hwe_blacks_X_chrom_final)
dim(snps_failed_hwe_blacks_X_chrom_final)
# [1] 1     1
colnames(snps_failed_hwe_blacks_X_chrom_final) <- 'SNP'
colnames(snps_failed_hwe_blacks_X_chrom_final) <- 'SNP'

# Add together all SNPs from autosomes and X chromosomes (females) that failed hwe
final_allethnicities <- rbind(snps_failed_hwe_mixed_autosomes_final,snps_failed_hwe_whites_autosomes_final,snps_failed_hwe_blacks_autosomes_final,
                              snps_failed_hwe_asians_autosomes_final,snps_failed_hwe_hispanics_autosomes_final, 
                              snps_failed_hwe_blacks_X_chrom_final)
dim(final_allethnicities)
# [1] 554     1

# save
# write.table(final_allethnicities, file = "snps_that_failed_hwe_in_all_ethnicities_in_autosomes_and_X_chrom_21dec2022.txt", quote = F, row.names = F, col.names = F)

# Plot of the heterozygosity rate distribution  -------------------------------------------------------
# Samples of white ancestry 

het_whites <- fread(paste0(dir,"whites_het_only_autosomes.het"), data.table = F)
het_whites$HET_RATE = (het_whites$"N(NM)" - het_whites$"O(HOM)")/het_whites$"N(NM)"

#pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/whites_heterozygosity_21dec2022.pdf")
hist(het_whites$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate (2000HIV of White ancestry)")
#dev.off()

# Plot of the heterozygosity rate as a function of the genotype missingness per sample
indmiss_whites_autosomes <- fread(paste0(dir,"whites_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes.imiss"), data.table = F)

head(indmiss_whites_autosomes)
head(het_whites )

indmiss_and_het_whites <- join(het_whites, indmiss_whites_autosomes, by = "IID")
head(indmiss_and_het_whites )

indmiss_and_het_whites_v2 <- indmiss_and_het_whites[,c("IID", "HET_RATE","F_MISS"),]
head(indmiss_and_het_whites_v2)

##calculate mean and standard deviation - heterozygosity rate
mean(indmiss_and_het_whites_v2$HET_RATE)
#[1] 0.1488218

sd(indmiss_and_het_whites_v2$HET_RATE)
#[1]  0.00157049
#3*sd
3* 0.00157049
#[1]  0.00471147

#mean+3sd
0.1488218+0.00471147
#[1]  0.1535333

#mean-3sd
0.1488218-0.00471147
#[1]   0.1441103

###############################################################################################################
##calculate mean and standard deviation - (proportion of missing SNPs per individual)
mean(indmiss_and_het_whites_v2$F_MISS)
#[1]   0.001032771

sd(indmiss_and_het_whites_v2$F_MISS)
#[1]0.001430021

#3*sd
3*0.001430021
#[1] 0.004290063

##mean+3sd
0.001032771+0.004290063
#[1] 0.005322834

##mean-3sd
0.001032771-0.004290063
#[1] -0.003257292

#####################################################################################

#Plot heterozygosity rate vs F MISS (proportion of missing SNPs per individual)
pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/whites_heterozygosity_vs_FMISS_21dec2022.pdf")
ggplot(indmiss_and_het_whites_v2, aes(x=F_MISS, y=HET_RATE)) +
  geom_point(stat="identity", fill="white")+
  labs(x="Proportion of missing SNPs per individual",y="Heterozygosity rate", size = 12)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=12),axis.title.y = element_text(size=11))+
  geom_hline(yintercept= 0.1488218 , color = "red")+
  geom_hline(yintercept= 0.1535333, color = "red", linetype="dashed")+
  geom_hline(yintercept=  0.1441103, color = "red", linetype="dashed")+
  geom_vline(xintercept= 0.001032771, color = "red")+
  scale_color_manual(values=c("red"))+
  ggtitle("2000HIV of White ancestry")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
het_whites$HET_RATE = (het_whites$"N(NM)" - het_whites$"O(HOM)")/het_whites$"N(NM)"
het_whites_fail = subset(het_whites, (het_whites$HET_RATE < mean(het_whites$HET_RATE)-3*sd(het_whites$HET_RATE)) | (het_whites$HET_RATE > mean(het_whites$HET_RATE)+3*sd(het_whites$HET_RATE)));
het_whites_fail$HET_DST = (het_whites_fail$HET_RATE-mean(het_whites$HET_RATE))/sd(het_whites$HET_RATE);
dim(het_whites_fail)
# [1] 32  8
# write.table(het_whites_fail, file = "whites_fail_het_qc_21dec2022.txt", row.names=FALSE, quote = F, sep = "\t")

# Plot of the heterozygosity rate distribution -------------------------------------------------------
# Samples of Black ancestry 

het_blacks <- fread(paste0(dir,"blacks_het_only_autosomes.het"), data.table = F)
het_blacks$HET_RATE = (het_blacks$"N(NM)" - het_blacks$"O(HOM)")/het_blacks$"N(NM)"

pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/blacks_heterozygosity_21dec2022.pdf")
hist(het_blacks$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate (2000HIV of Black ancestry)")
dev.off()

# Plot of the heterozygosity rate as a function of the genotype missingness per samples
indmiss_blacks_autosomes <- fread(paste0(dir,"blacks_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes.imiss"), data.table = F)

head(indmiss_blacks_autosomes)
head(het_blacks )

indmiss_and_het_blacks <- join(het_blacks, indmiss_blacks_autosomes, by = "IID")
head(indmiss_and_het_blacks )

indmiss_and_het_blacks_v2 <- indmiss_and_het_blacks[,c("IID", "HET_RATE","F_MISS"),]
head(indmiss_and_het_blacks_v2)

##calculate mean and standard deviation - heterozygosity rate
mean(indmiss_and_het_blacks_v2$HET_RATE)
#[1] 0.2091995

sd(indmiss_and_het_blacks_v2$HET_RATE)
#[1]  0.007312443
#3*sd
3* 0.007312443
#[1] 0.02193733

#mean+3sd
0.2091995+0.02193733
#[1]0.2311368

#mean-3sd
0.2091995-0.02193733
#[1]   0.1872622

###############################################################################################################
##calculate mean and standard deviation - (proportion of missing SNPs per individual)
mean(indmiss_and_het_blacks_v2$F_MISS)
#[1]0.001394114

sd(indmiss_and_het_blacks_v2$F_MISS)
#[1]0.001327276

#3*sd
3*0.001327276
#[1] 0.003981828

##mean+3sd
0.001394114+0.003981828
#[1] 0.005375942

##mean-3sd
0.001394114-0.003981828
#[1] -0.002587714

#####################################################################################

#Plot heterozygosity rate vs F MISS (proportion of missing SNPs per individual)
pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/blacks_heterozygosity_vs_FMISS_21dec2022.pdf")
ggplot(indmiss_and_het_blacks_v2, aes(x=F_MISS, y=HET_RATE)) +
  geom_point(stat="identity", fill="white")+
  labs(x="Proportion of missing SNPs per individual",y="Heterozygosity rate", size = 12)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=12),axis.title.y = element_text(size=11))+
  geom_hline(yintercept= 0.2091995, color = "red")+
  geom_hline(yintercept= 0.2311368, color = "red", linetype="dashed")+
  geom_hline(yintercept= 0.1872622, color = "red", linetype="dashed")+
  geom_vline(xintercept= 0.001394114, color = "red")+
  scale_color_manual(values=c("red"))+
  ggtitle("2000HIV of Black ancestry")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
het_blacks$HET_RATE = (het_blacks$"N(NM)" - het_blacks$"O(HOM)")/het_blacks$"N(NM)"
het_blacks_fail = subset(het_blacks, (het_blacks$HET_RATE < mean(het_blacks$HET_RATE)-3*sd(het_blacks$HET_RATE)) | (het_blacks$HET_RATE > mean(het_blacks$HET_RATE)+3*sd(het_blacks$HET_RATE)));
het_blacks_fail$HET_DST = (het_blacks_fail$HET_RATE-mean(het_blacks$HET_RATE))/sd(het_blacks$HET_RATE);
dim(het_blacks_fail)
# [1] 3  8
# write.table(het_blacks_fail, file = "blacks_fail_het_qc_21dec2022.txt", row.names=FALSE, quote = F, sep = "\t")

# Plot of the heterozygosity rate distribution -------------------------------------------------------
# Samples of Asian ancestry 

het_asians <- fread(paste0(dir,"asians_het_only_autosomes.het"), data.table = F)
het_asians$HET_RATE = (het_asians$"N(NM)" - het_asians$"O(HOM)")/het_asians$"N(NM)"

pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/asians_heterozygosity_21dec2022.pdf")
hist(het_asians$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate (2000HIV of Asian ancestry)")
dev.off()

# Plot of the heterozygosity rate as a function of the genotype missingness per samples
indmiss_asians_autosomes <- fread(paste0(dir,"asians_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes.imiss"), data.table = F)

head(indmiss_asians_autosomes)
head(het_asians )

indmiss_and_het_asians <- join(het_asians, indmiss_asians_autosomes, by = "IID")
head(indmiss_and_het_asians )

indmiss_and_het_asians_v2 <- indmiss_and_het_asians[,c("IID", "HET_RATE","F_MISS"),]
head(indmiss_and_het_asians_v2)

##calculate mean and standard deviation - heterozygosity rate
mean(indmiss_and_het_asians_v2$HET_RATE)
#[1]  0.1880095

sd(indmiss_and_het_asians_v2$HET_RATE)
#[1]  0.008604762
#3*sd
3* 0.008604762
#[1] 0.02581429

#mean+3sd
0.1880095+0.02581429
#[1]0.2138238

#mean-3sd
0.1880095-0.02581429
#[1]   0.1621952

###############################################################################################################
##calculate mean and standard deviation - (proportion of missing SNPs per individual)
mean(indmiss_and_het_asians_v2$F_MISS)
#[1]0.001252829

sd(indmiss_and_het_asians_v2$F_MISS)
#[1]0.001301487

#3*sd
3*0.001301487
#[1] 0.003904461

##mean+3sd
0.001252829+0.003904461
#[1]  0.00515729

##mean-3sd
0.001252829-0.003904461
#[1]  -0.002651632

#####################################################################################

#Plot heterozygosity rate vs F MISS (proportion of missing SNPs per individual)
pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/asians_heterozygosity_vs_FMISS_21dec2022.pdf")
ggplot(indmiss_and_het_asians_v2, aes(x=F_MISS, y=HET_RATE)) +
  geom_point(stat="identity", fill="white")+
  labs(x="Proportion of missing SNPs per individual",y="Heterozygosity rate", size = 12)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=12),axis.title.y = element_text(size=11))+
  geom_hline(yintercept= 0.1880095, color = "red")+
  geom_hline(yintercept= 0.2138238, color = "red", linetype="dashed")+
  geom_hline(yintercept= 0.1621952, color = "red", linetype="dashed")+
  geom_vline(xintercept= 0.001252829, color = "red")+
  scale_color_manual(values=c("red"))+
  ggtitle("2000HIV of Asian ancestry")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
het_asians$HET_RATE = (het_asians$"N(NM)" - het_asians$"O(HOM)")/het_asians$"N(NM)"
het_asians_fail = subset(het_asians, (het_asians$HET_RATE < mean(het_asians$HET_RATE)-3*sd(het_asians$HET_RATE)) | (het_asians$HET_RATE > mean(het_asians$HET_RATE)+3*sd(het_asians$HET_RATE)));
het_asians_fail$HET_DST = (het_asians_fail$HET_RATE-mean(het_asians$HET_RATE))/sd(het_asians$HET_RATE);
dim(het_asians_fail)
# [1] 1  8
# write.table(het_asians_fail, file = "asians_fail_het_qc_21dec2022.txt", row.names=FALSE, quote = F, sep = "\t")

# Plot of the heterozygosity rate distribution -------------------------------------------------------
# Samples of Hispanic ancestry 

het_hispanics <- fread(paste0(dir,"hispanics_het_only_autosomes.het"), data.table = F)
het_hispanics$HET_RATE = (het_hispanics$"N(NM)" - het_hispanics$"O(HOM)")/het_hispanics$"N(NM)"

pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/hispanics_heterozygosity_21dec2022.pdf")
hist(het_hispanics$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate (2000HIV of Hispanic ancestry)")
dev.off()

# Plot of the heterozygosity rate as a function of the genotype missingness per samples
indmiss_hispanics_autosomes <- fread(paste0(dir,"hispanics_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes.imiss"), data.table = F)

head(indmiss_hispanics_autosomes)
head(het_hispanics )

indmiss_and_het_hispanics <- join(het_hispanics, indmiss_hispanics_autosomes, by = "IID")
head(indmiss_and_het_hispanics )

indmiss_and_het_hispanics_v2 <- indmiss_and_het_hispanics[,c("IID", "HET_RATE","F_MISS"),]
head(indmiss_and_het_hispanics_v2)

##calculate mean and standard deviation - heterozygosity rate
mean(indmiss_and_het_hispanics_v2$HET_RATE)
#[1]   0.2182935

sd(indmiss_and_het_hispanics_v2$HET_RATE)
#[1] 0.008211279 
#3*sd
3* 0.008211279
#[1] 0.02463384

#mean+3sd
0.2182935+0.02463384
#[1] 0.2429273

#mean-3sd
0.2182935-0.02463384
#[1] 0.1936597

###############################################################################################################
##calculate mean and standard deviation - (proportion of missing SNPs per individual)
mean(indmiss_and_het_hispanics_v2$F_MISS)
#[1]0.001316813

sd(indmiss_and_het_hispanics_v2$F_MISS)
#[1] 0.00138652

#3*sd
3* 0.00138652
#[1] 0.00415956

##mean+3sd
0.001316813+0.00415956
#[1]0.005476373

##mean-3sd
0.001316813-0.00415956
#[1]  -0.002842747 

#####################################################################################

#Plot heterozygosity rate vs F MISS (proportion of missing SNPs per individual)
pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/hispanics_heterozygosity_vs_FMISS_21dec2022.pdf")
ggplot(indmiss_and_het_hispanics_v2, aes(x=F_MISS, y=HET_RATE)) +
  geom_point(stat="identity", fill="white")+
  labs(x="Proportion of missing SNPs per individual",y="Heterozygosity rate", size = 12)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=12),axis.title.y = element_text(size=11))+
  geom_hline(yintercept= 0.2182935, color = "red")+
  geom_hline(yintercept= 0.2429273, color = "red", linetype="dashed")+
  geom_hline(yintercept=  0.1936597, color = "red", linetype="dashed")+
  geom_vline(xintercept= 0.001316813, color = "red")+
  scale_color_manual(values=c("red"))+
  ggtitle("2000HIV of Hispanic ancestry")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
het_hispanics$HET_RATE = (het_hispanics$"N(NM)" - het_hispanics$"O(HOM)")/het_hispanics$"N(NM)"
het_hispanics_fail = subset(het_hispanics, (het_hispanics$HET_RATE < mean(het_hispanics$HET_RATE)-3*sd(het_hispanics$HET_RATE)) | (het_hispanics$HET_RATE > mean(het_hispanics$HET_RATE)+3*sd(het_hispanics$HET_RATE)));
het_hispanics_fail$HET_DST = (het_hispanics_fail$HET_RATE-mean(het_hispanics$HET_RATE))/sd(het_hispanics$HET_RATE);
dim(het_hispanics_fail)
# [1] 1  8
# write.table(het_hispanics_fail, file = "hispanics_fail_het_qc_21dec2022.txt", row.names=FALSE, quote = F, sep = "\t")

# Plot of the heterozygosity rate distribution -------------------------------------------------------
# Samples of native_american ancestry 

het_native_americans <- fread(paste0(dir,"native_americans_het_only_autosomes.het"), data.table = F)
het_native_americans$HET_RATE = (het_native_americans$"N(NM)" - het_native_americans$"O(HOM)")/het_native_americans$"N(NM)"

pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/native_americans_heterozygosity_21dec2022.pdf")
hist(het_native_americans$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate (2000HIV of Native american ancestry)")
dev.off()

# Plot of the heterozygosity rate as a function of the genotype missingness per samples
indmiss_native_americans_autosomes <- fread(paste0(dir,"native_americans_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes.imiss"), data.table = F)

head(indmiss_native_americans_autosomes)
head(het_native_americans )

indmiss_and_het_native_americans <- join(het_native_americans, indmiss_native_americans_autosomes, by = "IID")
head(indmiss_and_het_native_americans )

indmiss_and_het_native_americans_v2 <- indmiss_and_het_native_americans[,c("IID", "HET_RATE","F_MISS"),]
head(indmiss_and_het_native_americans_v2)

##calculate mean and standard deviation - heterozygosity rate
mean(indmiss_and_het_native_americans_v2$HET_RATE)
#[1]    0.4238872

sd(indmiss_and_het_native_americans_v2$HET_RATE)
#[1] 0.003963414 
#3*sd
3*0.003963414
#[1] 0.01189024

#mean+3sd
0.4238872+0.01189024
#[1] 0.4357774

#mean-3sd
0.4238872-0.01189024
#[1] 0.411997

###############################################################################################################
##calculate mean and standard deviation - (proportion of missing SNPs per individual)
mean(indmiss_and_het_native_americans_v2$F_MISS)
#[1] 0.0010563

sd(indmiss_and_het_native_americans_v2$F_MISS)
#[1] 0.0006524015

#3*sd
3*0.0006524015
#[1] 0.001957204

##mean+3sd
0.0010563+0.001957204
#[1] 0.003013504

##mean-3sd
0.0010563-0.001957204
#[1]-0.000900904 

#####################################################################################

#Plot heterozygosity rate vs F MISS (proportion of missing SNPs per individual)
pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/native_americans_heterozygosity_vs_FMISS_21dec2022.pdf")
ggplot(indmiss_and_het_native_americans_v2, aes(x=F_MISS, y=HET_RATE)) +
  geom_point(stat="identity", fill="white")+
  labs(x="Proportion of missing SNPs per individual",y="Heterozygosity rate", size = 12)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=12),axis.title.y = element_text(size=11))+
  geom_hline(yintercept= 0.4238872, color = "red")+
  geom_hline(yintercept= 0.4357774, color = "red", linetype="dashed")+
  geom_hline(yintercept=  0.411997, color = "red", linetype="dashed")+
  geom_vline(xintercept= 0.0010563, color = "red")+
  scale_color_manual(values=c("red"))+
  ggtitle("2000HIV of Native american ancestry")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
het_native_americans$HET_RATE = (het_native_americans$"N(NM)" - het_native_americans$"O(HOM)")/het_native_americans$"N(NM)"
het_native_americans_fail = subset(het_native_americans, (het_native_americans$HET_RATE < mean(het_native_americans$HET_RATE)-3*sd(het_native_americans$HET_RATE)) | (het_native_americans$HET_RATE > mean(het_native_americans$HET_RATE)+3*sd(het_native_americans$HET_RATE)));
het_native_americans_fail$HET_DST = (het_native_americans_fail$HET_RATE-mean(het_native_americans$HET_RATE))/sd(het_native_americans$HET_RATE);
dim(het_native_americans_fail)
# [1] 0  8
# write.table(het_native_americans_fail, file = "native_americans_fail_het_qc_21dec2022.txt", row.names=FALSE, quote = F, sep = "\t")

# Plot of the heterozygosity rate distribution -------------------------------------------------------
# Samples of mixed ancestry 

het_mixed <- fread(paste0(dir,"mixed_het_only_autosomes.het"), data.table = F)
het_mixed$HET_RATE = (het_mixed$"N(NM)" - het_mixed$"O(HOM)")/het_mixed$"N(NM)"

pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/mixed_heterozygosity_21dec2022.pdf")
hist(het_mixed$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate (2000HIV of Mixed ancestry)")
dev.off()

# Plot of the heterozygosity rate as a function of the genotype missingness per samples
indmiss_mixed_autosomes <- fread(paste0(dir,"mixed_mind0.025_hwe1minus6_geno0.05_updated_sex_new2_2000HIV_GSA796_only_autosomes.imiss"), data.table = F)

head(indmiss_mixed_autosomes)
head(het_mixed )

indmiss_and_het_mixed <- join(het_mixed, indmiss_mixed_autosomes, by = "IID")
head(indmiss_and_het_mixed )

indmiss_and_het_mixed_v2 <- indmiss_and_het_mixed[,c("IID", "HET_RATE","F_MISS"),]
head(indmiss_and_het_mixed_v2)

##calculate mean and standard deviation - heterozygosity rate
mean(indmiss_and_het_mixed_v2$HET_RATE)
#[1]     0.1729025

sd(indmiss_and_het_mixed_v2$HET_RATE)
#[1] 0.00582858 
#3*sd
3*0.00582858
#[1]  0.01748574

#mean+3sd
0.17290252+0.01748574
#[1] 0.1903883

#mean-3sd
0.17290252-0.01748574
#[1]  0.1554168

###############################################################################################################
##calculate mean and standard deviation - (proportion of missing SNPs per individual)
mean(indmiss_and_het_mixed_v2$F_MISS)
#[1] 0.001175659

sd(indmiss_and_het_mixed_v2$F_MISS)
#[1] 0.001253633

#3*sd
3*0.001253633
#[1] 0.003760899

##mean+3sd
0.001175659+0.003760899
#[1] 0.004936558

##mean-3sd
0.001175659-0.003760899
#[1] -0.00258524

#####################################################################################

#Plot heterozygosity rate vs F MISS (proportion of missing SNPs per individual)
pdf("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/plots/mixed_heterozygosity_vs_FMISS_21dec2022.pdf")
ggplot(indmiss_and_het_mixed_v2, aes(x=F_MISS, y=HET_RATE)) +
  geom_point(stat="identity", fill="white")+
  labs(x="Proportion of missing SNPs per individual",y="Heterozygosity rate", size = 12)+
  theme(panel.background = element_blank(),panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_text(size=12),axis.title.y = element_text(size=11))+
  geom_hline(yintercept= 0.1729025, color = "red")+
  geom_hline(yintercept= 0.1903883, color = "red", linetype="dashed")+
  geom_hline(yintercept=   0.1554168, color = "red", linetype="dashed")+
  geom_vline(xintercept= 0.001175659, color = "red")+
  scale_color_manual(values=c("red"))+
  ggtitle("2000HIV of Mixed ancestry")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# The following code generates a list of individuals who deviate more than 3 standard deviations from the heterozygosity rate mean.
het_mixed$HET_RATE = (het_mixed$"N(NM)" - het_mixed$"O(HOM)")/het_mixed$"N(NM)"
het_mixed_fail = subset(het_mixed, (het_mixed$HET_RATE < mean(het_mixed$HET_RATE)-3*sd(het_mixed$HET_RATE)) | (het_mixed$HET_RATE > mean(het_mixed$HET_RATE)+3*sd(het_mixed$HET_RATE)));
het_mixed_fail$HET_DST = (het_mixed_fail$HET_RATE-mean(het_mixed$HET_RATE))/sd(het_mixed$HET_RATE);
dim(het_mixed_fail)
# [1] 2  8
# write.table(het_mixed_fail, file = "mixed_fail_het_qc_21dec2022.txt", row.names=FALSE, quote = F, sep = "\t")

# Add all outliers due to excess heterozygosity in one file

het_whites_fail_final <- het_whites_fail[,1:2]
het_blacks_fail_final <- het_blacks_fail[,1:2]
het_asians_fail_final <- het_asians_fail[,1:2]
het_hispanics_fail_final <- het_hispanics_fail[,1:2]
het_mixed_fail_final <- het_mixed_fail[,1:2]

all_het_outliers <- rbind(het_whites_fail_final, het_blacks_fail_final,het_asians_fail_final,het_hispanics_fail_final,het_mixed_fail_final)
dim(all_het_outliers)
# [1] 39  2
head(all_het_outliers)
# write.table(all_het_outliers, file = "all_het_outliers.txt",row.names=FALSE, quote = F, sep = "\t")

# Additional plots ------------------------------------------------------------------------------------------------------------------------

# Check the distribution of HWE p-values of all SNPs.
#hwe_discovery <- fread(paste0(dir," "), data.table = F)
# pdf("f")
# hist(hwe_discovery[,9],main="Histogram HWE p-values of all SNPs\n from the 2000HIV Discovery cohort")
# dev.off()

# Selecting SNPs with HWE p-value below 0.00001, required for the plot below, allows to zoom in on strongly deviating SNPs. Do:
# awk '{if($9<0.00001) print$0}' file.hwe > output.hwe

# Generate plots to visualize the sex-check results.
#gender <- fread(paste0(dir,"myfile.sexcheck.csv"), data.table = F)

#pdf("output.pdf")
#hist(gender[,6],main="Gender", xlab="F")
#dev.off()

#pdf("output.pdf")
#male=subset(gender, gender$PEDSEX==1)
#hist(male[,6],main="Men",xlab="F")
#dev.off()

#pdf("output.pdf")
#female=subset(gender, gender$PEDSEX==2)
#hist(female[,6],main="Women",xlab="F")
#dev.off()

# THE END ------------------------------------------------------------------------------------------------------------------------------------
