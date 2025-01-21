# Readme
# Created on: 9 January 2023
# Source of raw 2000HIV genotype data: /mnt/workspace/1. Original data/1.Raw/2000HIV/2000HIV_GWAS_GSAData/GSA796/PLINK_220222_0252
# File names:
# GSA796.phenotype, GSA796.ped, GSA796.map, GSA796.bat, GSA796.script

# Project number: GSA796
# Samples included in analysis: 1917 
# Variants present in dataset: 725497 SNPs

# Genetic data from the 2000HIV cohort of all ethnicities were submitted for imputation on the TopMed server.
# - Number of samples submitted for imputation: 1864 
# - Number of SNPs submitted for imputation:: 493267 


# Number of samples of European ancestry (Whites), including family members: 1380

# Description - Filter and annotate the TopMed-imputed SNPs from the 2000 HIV cohort (all ethnicities)
# Only SNPs from chromosomes 1 - 22 -----------------------------------------------------------------------------------------------------------------------------------------

# Step 1 - Unzip the files using the password and move to folder '1_unzip'
inputFolder="/Volumes/diskAshur2/vicky/2000HIV_TopMed_imputed_data_jan2023/all_ethnicities/raw_vcf/"
midOut="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/1_unzip/"

for chr in {1..22}
do
mv $inputFolder/chr_${chr}/ $midOut/
done

mv $inputFolder/chr_X/ $midOut/

# Step 2 - Add proper header for BCFtools to work with ################################################################################################################################################
# write the header and then modify it using the following R script: michiganReheading.R

inputFolder="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/1_unzip/"
midOut="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/"

# Autosomes
for chr in {1..22}
do
mkdir -p $midOut/chr_${chr}/
done

# Chromosome X
mkdir -p $midOut/chr_X/

# Autosomes
midOut="${midOut}"
for chr in {1..22}
do
bcftools view -h ${inputFolder}chr_${chr}/chr${chr}.dose.vcf.gz > "${midOut}chr_${chr}/hdr.txt"
done

# Chromosome X
bcftools view -h ${inputFolder}chr_X/chrX.dose.vcf.gz > "${midOut}chr_X/hdr.txt"

# Autosomes 
for chr in {1..22}
do
midOut="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/chr_${chr}/"
cd $midOut

Rscript /Users/vickymatzaraki/Desktop/Bioinformatics/scripts_Linux/MICHIGAN_imputation_pipeline_vicky/6_filter_and_SNPannotation/michiganReheading.R  "${midOut}hdr.txt"
done

# Chromosome X
midOut="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/chr_X/"
cd $midOut

Rscript /Users/vickymatzaraki/Desktop/Bioinformatics/scripts_Linux/MICHIGAN_imputation_pipeline_vicky/6_filter_and_SNPannotation/michiganReheading.R  "${midOut}hdr.txt"

# Autosomes
midOut="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/chr_${chr}/"
output="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/"
for chr in {1..22}
do
vcfNewHeader="chr_${chr}_header.vcf.gz"
bcftools reheader ${inputFolder}chr_${chr}/chr${chr}.dose.vcf.gz \
--header "${midOut}hdr_michigan.txt" \
--output "${output}${vcfNewHeader}" 
done

output="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/"
for chr in {1..22}
do
tabix "${output}chr_${chr}_header.vcf.gz"
done

# Chromosome X
vcfNewHeader="chr_X_header.vcf.gz"
bcftools reheader ${inputFolder}chr_X/chrX.dose.vcf.gz \
--header "${midOut}hdr_michigan.txt" \
--output "${output}${vcfNewHeader}" 

tabix "${output}chr_X_header.vcf.gz"

# Step 3 - Subset and Filter per ethnicity group ################################################################################################################################################
# bcftools view -i'MAF[0]>=0.01 & (R2[0]>=0.3 || ER2[0]>=0.7)'
#### MAF >= 0.01 # For all ethnicites, do not filter based on MAF but later using the plink files per ethnicity group
#### R2 >= 0.3 # R2 is the squared correlation between input genotypes and imputed dosages. Rsq This is the estimated value of the squared correlation between imputed genotypes and true, unobserved genotypes. Since true genotypes are not available, this calculation is 
# based on the idea that poorly imputed genotype counts will shrink towards their expectations based on population allele frequencies alone;
#### ER2 >= 0.7 # With empirical r2 (ER2), the array genotypes are used as a ground truth to see how well the Markov model is doing at imputing genotypes. 
# The estimated r2 (R2) value uses allele frequency to calculate the errors in lieu of having genotypes to compare against. See: https://genome.sph.umich.edu/wiki/Minimac3_Info_File#EmpR.2C_EmpRsq

# Subset vcf files per ethnicity
# Autosomes
inputFolder="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/"
filterOut="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/3_whites/"

for chr in {1..22} X
do
mkdir -p $filterOut/chr_${chr}/
done

# Extract list of samples from vcf files using bcftools
# bcftools query -l /Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/2_header/chr_14_header.vcf.gz > sample.list.txt

# Whites (n = 1380 samples) -----------------------------------------------------------------------------
# Subset the vcf to only samples with white ancestry
# see: https://bioinformatics.stackexchange.com/questions/3477/how-to-subset-samples-from-a-vcf-file

# Autosomes and Chromosome X
for chr in {1..22} X
do
bcftools view -S updated_white_sample_list.txt ${inputFolder}chr_${chr}_header.vcf.gz  -Oz -o ${filterOut}chr_${chr}/chr_${chr}_2000HIV_whites_TopMed_imputed_all_ethnicities_GRCh38.vcf.gz  --force-samples
done

for chr in {1..22} X
do
tabix "${filterOut}chr_${chr}/chr_${chr}_2000HIV_whites_TopMed_imputed_all_ethnicities_GRCh38.vcf.gz"
done


# As with chromosomal sequences it is highly recommended (but not required) that the header include tags describing the contigs referred to in the VCF file.
# Nothing wrong, but rather an unexpected behaviour as:
# The contig tag is optional in VCFs
# Other bcftools filter commands don't need it
# I don't see any reason why subsetting on samples would particularly need it

# Do the filtering per ethnicity group ---------------------------------------------------------------------

# Whites (n = 1380) ----------------------------------------------------------------------------------------------------
mkdir -p /Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/4_filtered_whites/

inputFolder="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/3_whites/"
outPut="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/4_filtered_whites/"

for chr in {1..22} X
do
mkdir -p $outPut/chr_${chr}/
done


# Autosomes and Chromosome X
for chr in {1..22} X
do
bcftools view -i 'MAF[0]>=0.01 & (R2[0]>=0.3 || ER2[0]>=0.7)' ${inputFolder}chr_${chr}/chr_${chr}_2000HIV_whites_TopMed_imputed_all_ethnicities_GRCh38.vcf.gz \
--output-file "${outPut}chr_${chr}/chr_${chr}_2000HIV_whites_TopMed_imputed_filtered_GRCh38.vcf.gz" \
--output-type z 
done

for chr in {1..22} X
do
tabix "${outPut}chr_${chr}/chr_${chr}_2000HIV_whites_TopMed_imputed_filtered_GRCh38.vcf.gz"
done

# Number os SNPs on chrom X before filtering:  15331041
# Number os SNPs on chrom X after filtering: 360718



# Step 4 - Remove identifiers  ################################################################################################################################################
# Whites -----------------------------------------------------------------------------------------------------------

inputFolder="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/4_filtered_whites"
noIDout="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/5_noID_whites/"

# Autosomes and Chromosome X
for chr in {1..22} X
do
mkdir -p $noIDout/chr_${chr}/
done

for chr in {1..22} X
do
vcfFilteredNoID="chr_${chr}_2000HIV_whites_TopMed_imputed_filtered_GRCh38_noID.vcf.gz"
bcftools annotate ${inputFolder}/chr_${chr}/chr_${chr}_2000HIV_whites_TopMed_imputed_filtered_GRCh38.vcf.gz \
--output "${noIDout}chr_${chr}/${vcfFilteredNoID}" \
--output-type z \
--remove ID 
done

for chr in {1..22} X
do
tabix "${noIDout}/chr_${chr}/chr_${chr}_2000HIV_whites_TopMed_imputed_filtered_GRCh38_noID.vcf.gz"
done



# Step 5 -  VCF files: Change Chromosome Notation from chr to a numeric
# Source: https://www.biostars.org/p/98582/#332269


# whites -------------------------------------------------------------------------------------------------
inputFolder="/Volumes/diskAshur2/vicky/post_imputation_steps_all_ethnicities/"
updatedOut="${inputFolder}/6_updated_whites/"

for chr in {1..22} X
do
mkdir -p $updatedOut/chr_${chr}/
done

for i in {1..22} X
do
echo "chr$i $i" >> chr_name_conv.txt

vcfUpdated="updated_chr_${i}_2000HIV_whites_TopMed_imputed_filtered_GRCh38_noID.vcf.gz"
bcftools annotate \
--rename-chrs chr_name_conv.txt \
-Oz \
-o  "${updatedOut}chr_$i/${vcfUpdated}" \
${inputFolder}5_noID_whites/chr_${i}/chr_${i}_2000HIV_whites_TopMed_imputed_filtered_GRCh38_noID.vcf.gz  
done

for i in {1..22} X
do
tabix "${updatedOut}/chr_${i}/updated_chr_${i}_2000HIV_whites_TopMed_imputed_filtered_GRCh38_noID.vcf.gz"
done



# Step 6 - Annotate the SNPs using the TopMed panel 
# use the reference VCF file from NCBI which can be downloaded from their ftp, such as common SNPs, 
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz .
# However, if one wants all SNPs, see this, https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/All_20180418.vcf.gz .


# Whites -------------------------------------------------------------------------------------------------------------------------------
inputFolder="/Users/vickymatzaraki/Desktop/cluster_all_ethnicities_all_ethnicities/"
annotatedOut="${inputFolder}/7_annotated_whites/"

for chr in {1..22} X
do
mkdir -p $annotatedOut/chr_${chr}/
done

# Autosomes and Chromosome X
for chr in {1..22} X
do
vcfAnnotated="chr_${chr}_updated_2000HIV_whites_TopMed_imputed_filtered_GRCh38_annotated.vcf.gz"
bcftools annotate \
--annotations /Users/vickymatzaraki/Desktop/Bioinformatics/public_data/genome_references/human_9606_b151_GRCh38p7/All_20180418.vcf.gz  \
--columns ID \
--output "${annotatedOut}chr_${chr}/${vcfAnnotated}" \
--output-type z \
${inputFolder}6_updated_whites/chr_${chr}/updated_chr_${chr}_2000HIV_whites_TopMed_imputed_filtered_GRCh38_noID.vcf.gz 
done

# Create tabix
inputFolder="/Users/vickymatzaraki/Desktop/cluster_all_ethnicities_all_ethnicities/"
annotatedOut="${inputFolder}7_annotated_whites/"

for chr in {1..22} X
do
tabix "${annotatedOut}chr_${chr}/chr_${chr}_updated_2000HIV_whites_TopMed_imputed_filtered_GRCh38_annotated.vcf.gz"
done
