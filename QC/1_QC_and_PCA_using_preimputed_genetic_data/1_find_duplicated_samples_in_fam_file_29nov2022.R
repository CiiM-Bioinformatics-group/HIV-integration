# Readme
# Created on: 29 November 2022
# Source of raw genotype data: /mnt/workspace/1. Original data/1.Raw/2000HIV/2000HIV_GWAS_GSAData/GSA796/PLINK_220222_0252
# File names: GSA796.phenotype, GSA796.ped, GSA796.map, GSA796.bat, GSA796.script

# Load libraries ------------------------------------------------------------------------------------------------------------------------------------------------
library(data.table)
library(plyr)
library(ggrepel)
library(ggplot2)
library(dplyr)

# Aim 1: Find the duplicated FID and IID and replace to unique ids so we can distinguish them in the fam file -----------------------------------------------------------------

# Set working directory
setwd("Z:/2_Biostatistics_analyses/Vicky_Matzaraki/2000HIV_preimputed_genetic_data_all_ethnicities/raw_all_ethnicities/")
# Import the 2000HIV fam file
fam <- fread("2000HIV_GSA796_all_chromosomes.fam", data.table = F) 
dim(fam)
# [1] 1917    6
head(fam)

# Find duplicated IIDs
fam_dups <- subset(fam,duplicated(V2))
dim(fam_dups)
# [1] 3 6 # There are three duplicated values in the .fam file
# View(fam_dups)
#V1     V2  V3      V4 V5 V6
#1260 1134  OLV282  0  0  1 -9
#1453  982  OLV283  0  0  1 -9
#1772  802  OLV284  0  0  1 -9

# Convert the IID into unique rownames
rownames(fam) = make.names(fam$V1, unique=TRUE)
# View(fam)

## convert the rownames to columns
fam <- cbind(new_FID = rownames(fam), fam)

# Remove "X" from the FIDs
fam$new_FID <- gsub("X","",as.character(fam$new_FID))
head(fam)

# Check if there are duplicated IIDs
fam_dups2 <- subset(fam,duplicated(new_FID))
dim(fam_dups2)
# [1] 0 7

# Check the exact duplicated names for the following samples: OLV282,OLV283,OLV284
toMatch <- c("1134","982","802")  # 
matches <- unique(grep(paste(toMatch,collapse="|"), 
                       fam$new_FID, value=TRUE))
matches

# Select specific columns
fam2 <- fam[,c("new_FID","V2","V3","V4","V5","V6")]
head(fam2)

# save
# The duplicated values change from "1134" "982"  "802" to "1134.1" "982.1"  "802.1"
# write.table(fam2, file = "~/Desktop/2000HIV_GSA796_all_chromosomes.fam", quote = F, row.names = F, sep = "\t", col.names = F)

# The old and updated fam files are the same except the three duplicated sample ids: OLV282, OLV283, OLV284
