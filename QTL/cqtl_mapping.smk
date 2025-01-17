import yaml
import os
import json
import sys
import re
import string
import random
import traceback
import logging
import pandas as pd
import numpy as np

from glob import glob
from snakemake.logging import logger

wildcard_constraints:
    chr = "\d+",
    cov = "\w+"

# This snakefile is an adaptation of our full one from our pipeline "QTLmake".
# It is made specifically for our 2000HIV multi-omics integration cQTL project.
# The preimputation QC was done by Radboud University and is included in a separate folder.


# BEFORE RUNNING PIPELINE, EDIT THE LINE BELOW TO THE CORRECT 'tools' FOLDER
tools = "tools"
# should contain
# - TOPMed_freeze5
# - dbSNP
# - imputationbot
# - 1000G
####

# Based on cloned repository structure
scripts = workflow.basedir + "/workflow/scripts"
envs = workflow.basedir + "/workflow/envs"

## Parse config file
indir = config["i"] 								# path to input folder
outdir = config["o"] 								# path to output folder

cov_main = config.get("cov_main", "age,sex,bmi,season_sin,season_cos,lockdown,covid_vacc,covid,center_EMC,center_OLV")        # main effect covariates
cov_interact = config.get("cov_interact", "")       # interaction effect covariates
clump_kb = config.get("clump_kb", 500)              # clumping window size in kb for independent loci
clump_r2 = config.get("clump_r2", 0.05)             # clumping r2 threshold for independent loci
mapping_pval = config.get("mapping_pval", 1)        # p-value threshold for mapping
plot_manhattan = config.get("plot_manhattan", True) # whether to plot manhattan
plot_lambda = config.get("plot_lambda", False)       # whether to plot lambda
phenotype_type = config.get("phenotype", "cytokines") # what phenotype to use for mapping
outmap = outdir + "/" + phenotype_type              # output folder for mapping

# post imputation filtering thresholds
post_maf = config.get("post_maf", 0.01) # normally we use 0.05 but for 2000HIV we used 0.01 since we have a large cohort
post_hwe = config.get("post_hwe", 1e-12)
post_r2 = config.get("post_r2", 0.5)

shell("mkdir -p {outdir}")

if mode == "dosage":
    logger.info("DOSAGE MODE ACTIVATED - MAKE SURE THIS IS INTENTIONAL")
    CHR = glob_wildcards(indir+"/Genotype/dosage/chr{CHR}.txt").CHR
elif mode == "vcf":
    logger.info("VCF MODE ACTIVATED - MAKE SURE THIS IS INTENTIONAL")
    CHR = glob_wildcards(indir+"/Genotype/vcf/chr{CHR}.vcf.gz").CHR
else:
    CHR = range(1,23) 
#logger.info("Chromomosomes: %s " % CHR)


#Add the reference path for harmonizing
if refpanel == "hrc-r1.1":
	population = "eur"
	refpath = tools+"/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
	genome_build = "hg19"
elif refpanel == "topmed-r3":
	population  = "mixed"
	refpath = tools+"/TOPMed_freeze5/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz"
	genome_build = "hg38"

# Log and benchmark directories go into output folder
logdir = outdir+"/logs"
bmdir = outdir+"/benchmarks"

covariates_file = "%s/covariates.tsv" % (indir)
phenotype_path = "%s/%s/phenotype.tsv" % (indir, phenotype_type)

# Standard extensions for plink and plink2
ext = ["bed", "bim", "fam"]     # plink
ext2 = ["pgen", "psam", "pvar"] # plink2

# Start going through the input folder and printing some info
dirname = indir.split("/")[-1]
logger.info("\n> Input directory: %s\n" % indir)


COHORT = config["cohort"]
logger.info("> Cohort: %s\n" % COHORT)
stream = os.popen('zcat %s/Genotype/vcf/chr22.vcf.gz | grep -m1 "^#CHROM" | cut -f 10- | tr "\t" "\n"  | wc -l' % indir)

NUM_SAMPLES = int(stream.read().strip().split(" ")[0])
logger.info("> Num genotyped individuals: %s\n" % NUM_SAMPLES)

# Probably a hacky way to do this, maybe instead include them in a config file.
# Don't filter out cytokines anymore during the pipeline as it will then be
# expecting files that will never be created. Instead doing it pre-pipeline.
if os.path.isfile(phenotype_path):
    min_samples = 30 # Define the minimum number of repeats for filtering
    CYTOKINES = []
    df = pd.read_csv(phenotype_path, sep="\t", index_col=0)

    # exclude people that are not in the covariate file
    if os.path.isfile(covariates_file):
        stream = os.popen("""cat %s | head -1""" % covariates_file)
        cov_people = stream.read().strip().split("\t")[1:]
        existing_columns = [col for col in cov_people if col in df.columns]
        df = df[existing_columns]

    df_trim = df.copy()
    for cyt, values in df_trim.iterrows():
        if "rpmi" in cyt:
            df_trim = df_trim.drop(cyt)
            #logger.info("Cytokine %s excluded due to being RPMI" % cyt)
        elif "ctrl" in cyt:
            df_trim = df_trim.drop(cyt)
            #logger.info("Cytokine %s excluded due to being CTRL" % cyt)
        elif len(values[values.notnull()]) < min_samples:
            df_trim = df_trim.drop(cyt)
            logger.info("Cytokine %s excluded due to having less than %d samples" % (cyt, min_samples))
        else:
            CYTOKINES.append(cyt)
    
    logger.info(df_trim.head())

    shell("mkdir -p %s/phenotype" % outmap)
    df_trim.to_csv("%s/phenotype/filtered.tsv" % outmap, sep = "\t")
    logger.info("\n> Cytokines (%d): %s...\n" % (len(CYTOKINES), ", ".join(CYTOKINES[:10])))
    logger.info("When changing the cytokine input, please keep in mind to also delete cytokines_normalised.tsv if you actually want the cytokines file to be updated.\n")
else:
    logger.info("Cytokine file not present. This is fine is mapping is not being done.")
    CYTOKINES = "None"

if os.path.isfile(covariates_file):
    stream = os.popen("""cat %s | awk '{print $1}' | tail -n +2 """ % covariates_file)
    available_covariates = stream.read().strip().split("\n")
    logger.info("> Available covariates (%d): %s" % (len(available_covariates), ", ".join(available_covariates)))
    logger.info(" - Main covariates: " + ", ".join(cov_main.split(",")))
    diff_main = set(cov_main.split(",")) - set(available_covariates)
    if len(diff_main) > 1:
        logger.error("Main covariate %s not found in covariates file. Please check your covariates or config file.\n" % diff_main)
        sys.exit(1)

    if cov_interact != "":
        logger.info(" - Interact covariates: " + ", ".join(cov_interact.split(",")))
        diff_interact = set(cov_interact.split(",")) - set(available_covariates)
        if len(diff_interact) > 1:
            logger.error("Interact covariate %s not found in covariates file. Please check your covariates or config file." % diff_interact)
            sys.exit(1)
    else:
        logger.info(" - Chosen interact covariates: None")
    logger.info("\n")
else:
    logger.error("!!!! Covariates file not present !!!!")

# only keep "main" mapping if interaction is set to false
COVARIATES = ["main"]
if cov_interact != "":
    [COVARIATES.append(x) for x in cov_interact.split(",")]

onstart:
    print("##### Snakemake starting #####\n")
    print("\t Creating logs dir...\n")
    shell("mkdir -p {logdir}")
    print("\t Creating benchmarks dir...\n")
    shell("mkdir -p {bmdir}")

### START
# Rules that do not need to be submitted to cluster, when using --cluster (SLURM)
localrules: all, copy

def final_output():
    out = []

    out.append(expand(outmap+"/meta/main/{cyt}-.txt", cyt=CYTOKINES))
    out.append(expand(outmap+"/mapping/{cov}_1e-5.tsv", cov=COVARIATES))
    out.append(expand(outmap+"/mapping/{cov}_genomewide_loci.tsv", cov=COVARIATES))
    out.append(expand(outmap+"/mapping/{cov}_studywide_loci.tsv", cov=COVARIATES))

    # plots
    out.append(expand(outmap+"/plots/lead_snps/{cov}", cov=COVARIATES))
    out.append(expand(outmap+"/plots/locuszoom/{cov}", cov=COVARIATES))
    out.append(expand(outmap+"/plots/lead_snps/{cov}", cov=COVARIATES))

    if plot_manhattan:
        out.append(expand(outmap+"/plots/manhattan/{cov}/genomewide_manhattan.pdf", cov=COVARIATES))
    if plot_lambda:
        out.append(expand(outmap+"/plots/lambda/{cov}/lambda.png", cov=COVARIATES))
    return out

if False:
    logger.info("Final output files:")
    for x in final_output():
        logger.info("\n")
        if (type(x) == list):
            for y in x:
                logger.info(y)
        else:
            logger.info(x)

rule all:
    input:
        final_output()


rule copy:
    """ Soft link VCF, usually straight out of the imputation step """
    input:  	indir+"/Genotype/vcf/chr{chr}.vcf.gz"
    output: 	outdir+"/imputation/local/chr{chr}.dose.vcf.gz",
    log:    	logdir+"/copy/chr{chr}.log"
    benchmark: 	bmdir+"/copy/chr{chr}.txt"
    resources: 	runtime=10, mem_mb=500
    shell: "ln -sd {input} {output}"


rule filter:
    """ Filter VCF to only include SNPs with R2 >= 0.5 and MAF >= 0.05 and HWE >= 1e-12 """
    input:  	outdir+"/imputation/local/chr{chr}.dose.vcf.gz"  # either from copying vcf in vcf mode or output from imputation in non-vcf mode
    output: 	outdir+"/imputation/filter/chr{chr}.vcf.gz",
                #tbi = outdir+"/imputation/filter/chr{chr}.vcf.gz.tbi"
    params: 	maf = post_maf,
                hwe = post_hwe,
                r2 = post_r2,
                out = outdir+"/imputation/filter/chr{chr}"
    log:    	logdir+"/filter/chr{chr}.log"
    benchmark: 	bmdir+"/filter/chr{chr}.txt"
    resources: 	runtime=60, mem_mb=NUM_SAMPLES * 100  # needs to have unzipped vcf in memory which can be very big..
    shell: "plink2 --vcf {input} dosage=HDS --out {params.out} --export vcf vcf-dosage=HDS bgz \
            --extract-if-info 'R2 >= 0.5' --maf {params.maf} --hwe {params.hwe} \
            --rm-dup exclude-mismatch --set-all-var-ids \"@:#:\$r:\$a\" --new-id-max-allele-len 10000 missing"

rule index:
    """ Index VCF """
    input:  	rules.filter.output
    output: 	outdir+"/imputation/filter/chr{chr}.vcf.gz.csi"
    log:    	logdir+"/index/chr{chr}.log"
    benchmark: 	bmdir+"/index/chr{chr}.txt"
    conda:      envs+"/samtools.yaml"
    resources:  runtime=10, mem_mb=1000
    shell: "bcftools index {input}"

rule rsid:
    """ Add rsids to VCFs; final version of VCF before transitioning to other file types """
    input:  	vcf = rules.filter.output,
                csi = rules.index.output,
                dbsnp = tools+"/dbSNP/00-All.vcf.gz"
    output: 	outdir+"/imputation/rsid/chr{chr}.vcf.gz"
    log:    	logdir+"/rsid/chr{chr}.log"
    benchmark: 	bmdir+"/rsid/chr{chr}.txt"
    conda:      envs+"/samtools.yaml"
    resources:  runtime=360, mem_mb=1000
    shell: "bcftools annotate {input.vcf} -c =ID -a {input.dbsnp} -o {output}" # -l ID:unique

rule clumping_ref:
    """ Merge all chromosomes and create plink binary for clumping """
	input:		all = expand(rules.rsid.output, chr=CHR),
                first = expand(rules.rsid.output, chr=22)
	output: 	vcf = outdir+"/genotype/allchr.vcf",
                bed = outdir+"/genotype/allchr.bed"
	log:    	logdir+"/clumping_ref/log"
	benchmark: 	bmdir+"/clumping_ref/txt"
	resources:  runtime=120, mem_mb=60000
	shell: "zcat {input.first} | grep '^#' > {output.vcf} ;\
            zcat {input.all} | grep -h -v '^#' | awk 'BEGIN {{OFS=\"\t\"}} {{if (!($3 ~ /^chr/)) $3 = \"chr\" $3; print}}' >> {output.vcf} ;\
            plink2 --vcf {output.vcf} --make-bed --out {outdir}/genotype/allchr"

rule freq2:
    """ Get allele frequencies, needed for downstream MR analysis"""
    input:  	rules.rsid.output
    output: 	outdir+"/imputation/freq/chr{chr}.afreq"
    log:    	logdir+"/freq2/chr{chr}.log"
    benchmark: 	bmdir+"/freq2/chr{chr}.txt"
    conda:      envs+"/plink.yaml"
    resources: 	runtime=60, mem_mb=NUM_SAMPLES * 100
    shell: "plink2 --vcf {input} --freq --out {outdir}/imputation/freq/chr{wildcards.chr}"

rule alt_alleles: # doing it with bcftools means the vcf does not need to be gunzipped
    """ Get list of alt alleles """
    input:  	rules.rsid.output
    output: 	outdir+"/imputation/alt/chr{chr}.txt"
    log:    	logdir+"/alt_alleles/chr{chr}.log"
    benchmark: 	bmdir+"/alt_alleles/chr{chr}.txt"
    conda:      envs+"/samtools.yaml"
    resources:  runtime=5, mem_mb=1000
    shell: "bcftools query -f '%ID\t%ALT\n' {input} > {output}"

rule traw:
    """ Count alt alleles and make traw """
    input:      vcf = rules.rsid.output,
                alt_list = rules.alt_alleles.output
    output:     outdir+"/genotype/traw/chr{chr}.traw"
    log:        logdir+"/traw/chr{chr}.log"
    benchmark: 	bmdir+"/traw/chr{chr}.txt"
    resources:  runtime=5, mem_mb=lambda wildcards, input, attempt: (input.size//1000000)*(attempt*20)
    shell: "plink2 --vcf {input.vcf} dosage=HDS --export-allele {input.alt_list} \
            --export Av --out {outdir}/genotype/traw/chr{wildcards.chr}"

rule dosage:
    """ Format headers, add 'chr' to ID if missing, and exclude specific columns to get dosage format """
    input:  	rules.traw.output
    output: 	outdir+"/genotype/dosage/chr{chr}.txt"
    log:    	logdir+"/dosage/chr{chr}.log"
    benchmark: 	bmdir+"/dosage/chr{chr}.txt"
    resources:  runtime=5, mem_mb=lambda wildcards, input, attempt: (input.size//1000000)*(attempt*20)
    shell: "awk -F'\t' 'BEGIN {{OFS = FS}} NR==1 {{ for (i=1; i<=NF; i++) gsub(\"^([0-9]+_)+\", \"\", $i); $2 = \"SNP\" ;\
            print }}; NR>1{{if (!($2 ~ /^chr/)) $2 = \"chr\" $2; print }}' {input} | cut -f2,7- > {output}"

### GENOTYPE PROCESSING ENDS HERE ###
### PHENOTYPE PROCESSING STARTS HERE ###

rule normalise_phenotype:
    params:     inp = outmap+"/phenotype/filtered.tsv"  # since this file is updated every run, dont want to use it as input
    output:     outmap+"/phenotype/normalised.tsv"
    log:        logdir+"/normalise_phenotype/log"
    benchmark:  bmdir+"/normalise_phenotype.txt"
    resources:  runtime=max(int(len(CYTOKINES)/100), 1), mem_mb=5000
    shell: "Rscript {scripts}/normalise_phenotype.R {params.inp} {output} {outdir}/plots"

rule calc_studywide:
    """ Calculate study-wide significance threshold """
    input:      rules.normalise_phenotype.output #outmap+"/phenotype/normalised.tsv"
    output:     outmap+"/phenotype/multiple_testing.tsv"
    log:        logdir+"/calc_studywide/log"
    benchmark:  bmdir+"/calc_studywide.txt"
    resources:  runtime=10, mem_mb=1000
    shell: "Rscript {scripts}/calc_studywide.R {input} {output}"

### PHENOTYPE PROCESSING ENDS HERE ###
### MAPPING STARTS HERE ###

rule mapping:
    """ Run mapping; for non-interaction model and for interaction models using each covariate as interaction term """
    input:  	snp = rules.dosage.output,
                cyt = rules.normalise_phenotype.output
    output: 	expand(outmap+"/mapping-chr/{{cov}}/chr{{chr}}/{cyt}.tsv", cyt=CYTOKINES)
    params: 	cov = covariates_file,
                out = outmap+"/mapping-chr/{cov}/chr{chr}",
                p_val = mapping_pval
    log:    	logdir+"/mapping/chr{chr}/{cov}.log"
    benchmark: 	bmdir+"/mapping/chr{chr}/{cov}.txt"
    resources:  runtime=int(len(CYTOKINES)/800 * (6*NUM_SAMPLES)), mem_mb=lambda wildcards, input, attempt: (input.size//1000000)*(attempt*50)
    shell: "Rscript {scripts}/run_matrixEQTL.R {input.snp} {params.cov} {input.cyt} \
            {params.out} {params.p_val} {cov_main} {wildcards.cov} cytokines" # {cov_interact}

rule calc_lambda:
    input:		expand(outmap+"/mapping-chr/{{cov}}/chr{chr}/{{cyt}}.tsv", chr=CHR)
    output:     outmap+"/plots/lambda/{cov}/val/{cyt}.txt"
    log:        logdir+"/calc_lambda/{cov}/{cyt}.log"
    benchmark: 	bmdir+"/calc_lambda/{cov}/{cyt}.txt"
    resources: 	runtime=15, mem_mb=5000
    shell: "Rscript {scripts}/calc_lambda.R {output} {input}"  # keep input as last argument

rule plot_lambda:
    input:		expand(outmap+"/plots/lambda/{{cov}}/val/{cyt}.txt", cyt=CYTOKINES),
    output:     outmap+"/plots/lambda/{cov}/lambda.png"
    params:     mapping = outmap+"/mapping-chr/{cov}"
    log:        logdir+"/plot_lambda/{cov}.log"
    benchmark: 	bmdir+"/plot_lambda/{cov}.txt"
    resources: 	runtime=10, mem_mb=3000
    shell: "Rscript {scripts}/plot_lambda.R {output} {params.mapping} {input}"  # keep input as last argument

rule merge_chromosomes:
    """ Merge all results into one file per chromosome """
    input:  	expand(outmap+"/mapping-chr/{{cov}}/chr{{chr}}/{cyt}.tsv", cyt=CYTOKINES)
    output: 	outmap+"/mapping/{cov}/chr{chr}.tsv"
    log:    	logdir+"/merge_chromosomes/{cov}/chr{chr}.log"
    benchmark: 	bmdir+"/merge_chromosomes/{cov}/chr{chr}.txt"
    resources:  runtime=15, mem_mb=1000
    shell: "awk 'NR == 1 || FNR > 1' {input} > {output}" # only takes header from first file

rule get_significant:
    """ Extract only nominal significant results, makes making manhattan plot faster by doing it like this first"""
	input:		rules.merge_chromosomes.output
	output: 	outmap+"/mapping/{cov}/nominal/chr{chr}.tsv"
	log:    	logdir+"/get_significant/{cov}/chr{chr}.log"
	benchmark: 	bmdir+"/get_significant/{cov}/chr{chr}.txt"
	resources:  runtime=lambda wildcards, input, attempt: NUM_SAMPLES * len(CYTOKINES) / 20000 + (attempt*60), mem_mb=20000
	shell: "awk 'NR == 1 || $5 < 0.05' {input} > {output}"

rule merge_significant:
    """ Merge all significant results"""
	input:		expand(outmap+"/mapping/{{cov}}/nominal/chr{chr}.tsv", chr=CHR)
	output: 	outmap+"/mapping/{cov}_nominal.tsv"
	log:    	logdir+"/merge_significant/{cov}.log"
	benchmark: 	bmdir+"/merge_significant/{cov}.txt"
	resources:  runtime=lambda wildcards, input, attempt: NUM_SAMPLES * len(CYTOKINES) / 40000 + (attempt*60), mem_mb=30000
	shell: "awk 'NR == 1 || FNR > 1' {input} > {output}" # only takes header from first file

rule get_1e5_significant:
    """ Merge all cytokines and extract only 1e-5 significant results"""
	input: 		rules.merge_significant.output
	output: 	outmap+"/mapping/{cov}_1e-5.tsv"
	log:    	logdir+"/get_significant/{cov}.log"
	benchmark: 	bmdir+"/get_significant/{cov}.txt"
	resources:  runtime=60, mem_mb=6000
	shell: "awk 'NR==1 || $5 < 1e-5' {input} > {output}"

rule get_genomewide_significant:
    """ Merge all cytokines and extract only genome-wide significant results"""
	input: 		rules.merge_significant.output
	output: 	outmap+"/mapping/{cov}_genomewide.tsv"
	log:    	logdir+"/get_significant/{cov}.log"
	benchmark: 	bmdir+"/get_significant/{cov}.txt"
	resources:  runtime=60, mem_mb=6000
	shell: "awk 'NR==1 || $5 < 5e-8' {input} > {output}"

rule get_studywide:
    """ Extract only study-wide significant results"""
    input:		inp = outmap+"/mapping/{cov}_genomewide.tsv",
                threshold = rules.calc_studywide.output
    output: 	outmap+"/mapping/{cov}_studywide.tsv"
    log:    	logdir+"/get_studywide/{cov}.log"
    benchmark: 	bmdir+"/get_studywide/{cov}.txt"
    resources:  runtime=lambda wildcards, input, attempt: NUM_SAMPLES * len(CYTOKINES) / 60000 + (attempt*60), mem_mb=6000
    shell: "threshold=$(tail -1 {input.threshold} | awk '{{print $3}}'); awk -v threshold=\"$threshold\" 'NR == 1 || $5 < threshold' {input.inp} > {output}"

rule plot_manhattan:
    """ Plot manhattan plot """
	input: 		qtl = expand(outmap+"/mapping/{{cov}}/nominal/chr{chr}.tsv", chr=CHR)
	output: 	outmap+"/plots/manhattan/{cov}/genomewide_manhattan.pdf"
	params: 	out = outmap+"/plots/manhattan/{cov}"
	log:    	logdir+"/plot_manhattan/{cov}.log"
	benchmark: 	bmdir+"/plot_manhattan/{cov}.txt"
	resources:  runtime=120, mem_mb=lambda wildcards, input, attempt: max(200000, input.size//(700000/attempt))
	shell: "Rscript {scripts}/plot_manhattan.R {params.out} {wildcards.cov} {genome_build} FALSE {output} 5e-8 {input.qtl}"

rule table_loci:
    """ Get loci from all cytokines and make table """
    input:      gw = rules.get_genomewide_significant.output, # and mapping output
                bed = rules.clumping_ref.output.bed
    output:     outmap+"/mapping/{cov}_genomewide_loci.tsv"
    log:        logdir+"/table_loci/{cov}.log"
    benchmark:  bmdir+"/table_loci/{cov}.txt"
    resources:  runtime=60, mem_mb=3000
    shell: "Rscript {scripts}/table_loci.R {input.gw} {clump_kb} {clump_r2} {output} \
            {outmap}/plots {wildcards.cov} {outmap}/mapping-chr {input.bed} {outmap}/mapping-chr/main/chr1/mapped_samples.txt"

rule table_loci_studywide:
    """ Get loci from all cytokines and make table """
    input:      gw = rules.get_studywide.output, # and mapping output
                bed = rules.clumping_ref.output.bed
    output:     outmap+"/mapping/{cov}_studywide_loci.tsv"
    log:        logdir+"/table_loci/{cov}.log"
    benchmark:  bmdir+"/table_loci/{cov}.txt"
    resources:  runtime=60, mem_mb=3000
    shell: "Rscript {scripts}/table_loci.R {input.gw} {clump_kb} {clump_r2} {output} \
            {outmap}/plots {wildcards.cov} {outmap}/mapping-chr {input.bed} {outmap}/mapping-chr/main/chr1/mapped_samples.txt"
