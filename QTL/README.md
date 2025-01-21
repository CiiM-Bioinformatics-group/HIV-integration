# QTL mapping


The snakefile contains all the steps in chronological order. Firstly, some post-imputation genotype processing (outcome of the QC steps), such as creating dosage files. Then phenotype processing, such as normalisation, as well as calculating the study-wide P-value. Then finally the QTL mapping, followed by clumping to find individual loci, and plotting the results.