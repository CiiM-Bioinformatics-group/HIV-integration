# Mendelian Randomisation


The snakefile contains all the steps in chronological order for the MR analysis, e.g. getting the instrumental variables from the exposure, validating them, getting the same SNPs from the outcome, performing MR, and lastly performing sensitivity checks. At the same time, all these steps are being done in the opposite direction, in order to then exclude any significant MR results that also turned out to be significant in the opposite direction. 