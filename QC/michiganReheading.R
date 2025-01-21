#!/usr/bin/Rscript
#### 20 Sep 2016
####
#### Adding information to the header of michigan imputed vcf files.
## michiganReheading.R

args <- commandArgs(TRUE)
headerFile <- args[1]
outputFileName <- "hdr_michigan.txt"
#headerFile <- "hdr.txt"
head <- readLines(headerFile)
#thead <- readLines(outputFileName)

genoytpeRows <- c("##FILTER=<ID=GENOTYPED,Description=\"Site was genotyped\">", 
					"##FILTER=<ID=GENOTYPED_ONLY,Description=\"Site was genotyped only\">")


### search for the position
if(sum(!grepl(head, pattern="##FILTER=<ID=GENOTYPED,Description=\"Site was genotyped\">")) == length(head)) {
	upper <- head[1:grep(pattern="##INFO=<ID=ER2*", head)]

	bottom <- head[(grep(pattern="##INFO=<ID=ER2*", head)+1):length(head)]

	completeHeader <- c(upper, genoytpeRows, bottom)
	cat(completeHeader, file=outputFileName, sep="\n")
	
	} else {

	cat(head, file=outputFileName, sep="\n")
}
cat("New header out. hdr_michigan.txt")
### include 
##FILTER=<ID=GENOTYPED,Description="Site was genotyped">
##FILTER=<ID=GENOTYPED_ONLY,Description="Site was genotyped only">

### save new header and use bcftools to re-header the vcf file.  
