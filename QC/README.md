# Genotype preprocessing, QC and imputation


"1_QC_and_PCA_using_preimputed_genetic_data" contains the QC steps of the pre-imputed dataset, and the PCA analysis using the 1000Genomes data, as well as two additional R scripts with small things like preparing plots of text files for the QC steps. 

"2_post_imputation_QC_steps.sh" takes the imputed VCFs, extracts the European individuals, and performs some more QC steps. michiganReheading.R accompanies this script.

The outcome of the 2_post_imputation_QC_steps.sh is what is used as input for the QTL snakefile found in the QTL folder (../QTL).