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

# Based on cloned repository structure
scripts = workflow.basedir + "/workflow/scripts/mr"
#envs = workflow.basedir + "/workflow/envs"

#### start parameters ####
rootdir = config.get("dir", ".")
cohort = config.get("cohort", "2000HIV-EU-discovery")
cov = config.get("cov", "main")
exposure = config.get("exposure", "proteins") # metabolites, expression
outcome = config.get("outcome", "cytokines")
pval = config.get("pval", "1e-5")

do_validation = config.get("do_validation", "False")
val_pval_exposure = 0.05 #config.get("val_pval_exposure", "0.05")
val_pval_outcome = 0.05 #config.get("val_pval_outcome", "0.1")
val_cohort = config.get("val_cohort", "2000HIV-EU-validation")
#### end parameters ####

outdir = rootdir+"/"+cohort+"/mr/"+cov
shell("mkdir -p {outdir}")

# Log and benchmark directories go into output folder
logdir = outdir+"/logs/"+exposure
bmdir = outdir+"/benchmarks/"+exposure

# Parse phenotype files
exposure_path = rootdir+"/"+cohort+"/"+exposure+"/phenotype/normalised.tsv"
if os.path.isfile(exposure_path):
    EXPOSURE = []
    df = pd.read_csv(exposure_path, sep="\t", index_col=0, header=0)
    for pheno, values in df.iterrows():
        # <<< filtering here if needed >>>
        EXPOSURE.append(str(pheno))
    logger.info("\nExposure -> %s (%d): %s...\n" % (exposure.capitalize(), len(EXPOSURE), ", ".join(EXPOSURE[:10]))) 
else:
    logger.info("Exposure (%s) file not present." % exposure.capitalize())
    stop()

outcome_path = rootdir+"/"+cohort+"/"+outcome+"/phenotype/normalised.tsv"
if os.path.isfile(outcome_path):
    OUTCOME = []
    df = pd.read_csv(outcome_path, sep="\t", index_col=0, header=0)
    for pheno, values in df.iterrows():
        # <<< filtering here if needed >>>
        if outcome == "cytokines" and "rpmi" in pheno:
            continue
        OUTCOME.append(str(pheno))
    logger.info("Outcome -> %s (%d): %s...\n" % (outcome.capitalize(), len(OUTCOME), ", ".join(OUTCOME[:10]))) 
else:
    logger.info("Outcome (%s) file not present." % outcome.capitalize())
    stop()

CHR = range(1, 22)

# testing
CHR = 21
#EXPOSURE = EXPOSURE[:10]
#EXPOSURE = ["C9H9N", "H2O3S2", "H2O3S", "H2O4S"]

onstart:
    print("##### Snakemake starting #####\n")
    print("\t Creating logs dir...\n")
    shell("mkdir -p {logdir}")
    print("\t Creating benchmarks dir...\n")
    shell("mkdir -p {bmdir}")

### START OF PIPELINE ###
localrules: all

rule all:
    input:  
        expand(outdir+"/"+outcome+"_on_"+exposure+"/per_outcome/{e}.tsv", e = EXPOSURE),
        expand(outdir+"/"+exposure+"_on_"+outcome+"/per_outcome/{o}.tsv", o = OUTCOME),
        outdir+"/"+exposure+"_on_"+outcome+"/significant_mr_checked_pval.tsv",
        outdir+"/"+outcome+"_on_"+exposure+"/significant_mr_checked_pval.tsv"
    
rule calculate_pval:
    input:      one = rootdir+"/"+cohort+"/"+exposure+"/phenotype/multiple_testing.tsv",
                two = rootdir+"/"+cohort+"/"+outcome+"/phenotype/multiple_testing.tsv"
    output:     outdir+"/"+exposure+"_on_"+outcome+"/pvalue.txt"
    log:        logdir+"/calculate_pval/log"
    benchmark:  bmdir+"/calculate_pval.txt"
    resources:  mem_mb = 10000, runtime = 60
    shell: "Rscript {scripts}/calculate_pval.R {input.one} {input.two} {output}"

rule clump_exposure:
    input:      rootdir+"/"+cohort+"/"+exposure+"/mapping/main_"+pval+".tsv",
                expand(rootdir+"/"+cohort+"/imputation/freq/chr{chr}.afreq", chr = range(1, 22))
    output:     expand(outdir+"/clump/"+exposure+"/{e}.tsv", e = EXPOSURE)
    params:     ref = rootdir+"/"+cohort+"/genotype/allchr",
                pval_val = val_pval_exposure
    log:        logdir+"/clump_exposure/log"
    benchmark:  bmdir+"/clump_exposure.txt"
    resources:  mem_mb = int(len(EXPOSURE)*20)+10000, runtime = max(60, int(len(EXPOSURE)/2)+20), 
    shell: "Rscript {scripts}/genomewide_clumping.R {rootdir} {cohort} {cov} {exposure} {params.ref} {pval} {outdir}/clump/{exposure} {do_validation} {val_cohort} {params.pval_val}"

rule clump_outcome:
    input:      rootdir+"/"+cohort+"/"+outcome+"/mapping/main_"+pval+".tsv",
                expand(rootdir+"/"+cohort+"/imputation/freq/chr{chr}.afreq", chr = range(1, 22))
    output:     expand(outdir+"/clump/"+outcome+"/{o}.tsv", o = OUTCOME)
    params:     ref = rootdir+"/"+cohort+"/genotype/allchr",
                pval_val = val_pval_outcome
    log:        logdir+"/clump_outcome/log"
    benchmark:  bmdir+"/clump_outcome.txt"
    resources:  mem_mb = int(len(OUTCOME)*20)+10000, runtime = max(60, int(len(OUTCOME)/2)+20)
    shell: "Rscript {scripts}/genomewide_clumping.R {rootdir} {cohort} {cov} {outcome} {params.ref} {pval} {outdir}/clump/{outcome} {do_validation} {val_cohort} {params.pval_val}"

rule get_snps:
    # get genetic IVs from exposure from outcome
    input:      outdir+"/clump/"+exposure+"/{e}.tsv"
    output:     outdir+"/"+exposure+"_on_"+outcome+"/data/{e}_snps.tsv"
    params:     exposure = exposure,
                outcome = outcome,
                indir = outdir+"/"+exposure+"_on_"+outcome+"/data"
    log:        logdir+"/get_snps/{e}.log"
    benchmark:  bmdir+"/get_snps/{e}.txt"
    resources:  mem_mb = 10000, runtime = int(len(OUTCOME)/2)+30
    shell: "Rscript {scripts}/get_snps.R {input} {params.exposure} {params.outcome} {params.indir} {rootdir} {cohort} {cov}"

use rule get_snps as get_snps_bidirectional with:
    # get genetic IVs from exposure from outcome
    input:      outdir+"/clump/"+outcome+"/{o}.tsv"
    output:     outdir+"/"+outcome+"_on_"+exposure+"/data/{o}_snps.tsv"
    params:     exposure = outcome,
                outcome = exposure,
                indir = outdir+"/"+outcome+"_on_"+exposure+"/data"
    log:        logdir+"/get_snps_bidirectional/{o}.log"
    benchmark:  bmdir+"/get_snps_bidirectional/{o}.txt"
    resources:  mem_mb = 10000, runtime = int(len(EXPOSURE)/2)+30

# TODO, fix parameters inside script
# TODO: move check for pleiotropy to after clumping
rule mr:
    input:      file1 = outdir+"/clump/"+exposure+"/{e}.tsv",  # exposure snps are genetic instruments,
                file2 = outdir+"/"+exposure+"_on_"+outcome+"/data/{e}_snps.tsv"
    output:     outdir+"/"+exposure+"_on_"+outcome+"/per_exposure/{e}/all_mr.tsv"
    params:     exposure = exposure,
                outcome = outcome,
                outpath = outdir+"/"+exposure+"_on_"+outcome+"/per_exposure/{e}"
    log:        logdir+"/mr/{e}.log"
    benchmark:  bmdir+"/mr/{e}.txt"
    resources:  mem_mb = 10000, runtime = int(len(OUTCOME)/20)+30#, qos = "long"
    shell: "Rscript {scripts}/mr_indiv.R {input.file1} {input.file2} {params.exposure} {params.outcome} {params.outpath}"

use rule mr as mr_bidirectional with:   
    input:      file1 = outdir+"/clump/"+outcome+"/{o}.tsv",   # now outcome snps are genetic instruments,
                file2 = outdir+"/"+outcome+"_on_"+exposure+"/data/{o}_snps.tsv"
    output:     outdir+"/"+outcome+"_on_"+exposure+"/per_exposure/{o}/all_mr.tsv"
    params:     exposure = outcome,
                outcome = exposure,
                outpath = outdir+"/"+outcome+"_on_"+exposure+"/per_exposure/{o}"
    log:        logdir+"/mr_bidirectional/{o}.log"
    benchmark:  bmdir+"/mr_bidirectional/{o}.txt"
    resources:  mem_mb = 10000, runtime = int(len(EXPOSURE)/20)+30#, qos = "long"

rule combine_mr:
    input:      expand(rules.mr.output, e = EXPOSURE)
    output:     outdir+"/"+exposure+"_on_"+outcome+"/significant_mr.tsv"
    params:     indir = outdir+"/"+exposure+"_on_"+outcome+"/per_exposure",
                outdir = outdir+"/"+exposure+"_on_"+outcome,
                pheno1 = exposure,
                pheno2 = outcome
    log:        logdir+"/combine_mr/log"
    benchmark:  bmdir+"/combine_mr.txt"
    resources:  mem_mb = 10000, runtime = int(len(EXPOSURE)/20+30)
    shell: "Rscript {scripts}/combine_mr.R {params.indir} {params.outdir} {params.pheno1} {params.pheno2}"

use rule combine_mr as combine_mr_bidirectional with:
    input:      expand(rules.mr_bidirectional.output, o = OUTCOME)
    output:     outdir+"/"+outcome+"_on_"+exposure+"/significant_mr.tsv"
    params:     indir = outdir+"/"+outcome+"_on_"+exposure+"/per_exposure",
                outdir = outdir+"/"+outcome+"_on_"+exposure,
                pheno1 = outcome,
                pheno2 = exposure
    log:        logdir+"/combine_mr_bidirectional/log"
    benchmark:  bmdir+"/combine_mr_bidirectional.txt"
    resources:  mem_mb = 10000, runtime = int(len(OUTCOME)/20+30)

# check if mr is bidirectional, exclude if so
rule check_direction:
    input:      one = outdir+"/"+exposure+"_on_"+outcome+"/significant_mr.tsv",
                two = outdir+"/"+outcome+"_on_"+exposure+"/significant_mr.tsv"
    output:     outdir+"/"+exposure+"_on_"+outcome+"/significant_mr_checked.tsv"
    log:        logdir+"/check_direction/log"
    benchmark:  bmdir+"/check_direction.txt"
    resources:  mem_mb = 10000, runtime = 60
    shell: "Rscript {scripts}/check_direction.R {output} {input.one} {input.two}"

use rule check_direction as check_direction_bidirectional with:
    input:      one = outdir+"/"+outcome+"_on_"+exposure+"/significant_mr.tsv",
                two = outdir+"/"+exposure+"_on_"+outcome+"/significant_mr.tsv"
    output:     outdir+"/"+outcome+"_on_"+exposure+"/significant_mr_checked.tsv"
    log:        logdir+"/check_direction_bidirectional/log"
    benchmark:  bmdir+"/check_direction_bidirectional.txt"

rule check_pval:
    input:      mr = outdir+"/"+exposure+"_on_"+outcome+"/significant_mr_checked.tsv",
                pval = rules.calculate_pval.output
    output:     outdir+"/"+exposure+"_on_"+outcome+"/significant_mr_checked_pval.tsv"
    log:        logdir+"/check_pval/log"
    benchmark:  bmdir+"/check_pval.txt"
    resources:  mem_mb = 10000, runtime = 60
    shell: "Rscript {scripts}/check_pval.R {output} {input.mr} {input.pval}"

use rule check_pval as check_pval_bidirectional with:
    input:      mr = outdir+"/"+outcome+"_on_"+exposure+"/significant_mr_checked.tsv",
                pval = rules.calculate_pval.output
    output:     outdir+"/"+outcome+"_on_"+exposure+"/significant_mr_checked_pval.tsv"
    log:        logdir+"/check_pval_bidirectional/log"
    benchmark:  bmdir+"/check_pval_bidirectional.txt"

# combine per outcome
rule split_outcome:
    input:      rules.check_direction.output
    output:     expand(outdir+"/"+exposure+"_on_"+outcome+"/per_outcome/{o}.tsv", o = OUTCOME)
    params:     outdir = outdir+"/"+exposure+"_on_"+outcome+"/per_outcome",
                outcome = outcome
    log:        logdir+"/combine_mr_outcome/log"
    benchmark:  bmdir+"/combine_mr_outcome.txt"
    resources:  mem_mb = 10000, runtime = int(len(EXPOSURE)/60+30)
    shell: "Rscript {scripts}/split_outcome.R {input} {params.outdir} {params.outcome}"

use rule split_outcome as split_outcome_bidirectional with:
    input:      rules.check_direction_bidirectional.output
    output:     expand(outdir+"/"+outcome+"_on_"+exposure+"/per_outcome/{e}.tsv", e = EXPOSURE)
    params:     outdir = outdir+"/"+outcome+"_on_"+exposure+"/per_outcome",
                outcome = exposure
    log:        logdir+"/combine_mr_outcome_bidirectional/log"
    benchmark:  bmdir+"/combine_mr_outcome_bidirectional.txt"
    resources:  mem_mb = 10000, runtime = int(len(OUTCOME)/60+30)