# Snakemake file for ChIP-Seq PE analysis
# The pipeline runs with the command "snakemake --use-conda [--core] ", within the folder containing the Snakefile, after activation of the global_env.yaml

###############
# Libraries
###############

import os
import pandas as pd
from snakemake.utils import validate, min_version

###########
# Container
###########
singularity: "shub://truatpasteurdotfr/singularity-docker-miniconda"

#############################################
# Configuration and sample sheets
#############################################

configfile: "config.yaml"

WORKING_DIR         = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR          = config["result_dir"]      # what you want to keep

GENOME_FASTA_URL    = config["refs"]["genome_url"]
GENOME_FASTA_FILE   = os.path.basename(config["refs"]["genome_url"])
GFF_URL             = config["refs"]["gff_url"]
GFF_FILE            = os.path.basename(config["refs"]["gff_url"])

TOTALCORES          = 16                             #check this via 'grep -c processor /proc/cpuinfo'

###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_samples_per_treatment(input_df="units.tsv",colsamples="sample",coltreatment="condition",treatment="control"):
    """This function returns a list of samples that correspond to the same experimental condition"""
    df = pd.read_table(input_df)
    df = df.loc[df[coltreatment] == treatment]
    filtered_samples = df[colsamples].tolist()
    return filtered_samples

def is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv, it is used by the function get_trimmed_reads to define the samples used by the align rule"""
    return pd.isnull(units.loc[(sample), "fq2"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of a given sample """
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(WORKING_DIR + "trimmed/{sample}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return WORKING_DIR + "trimmed/{sample}.fastq.gz".format(**wildcards)

##############
# Samples and conditions
##############

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

CASES = get_samples_per_treatment(treatment="treatment")
CONTROLS = get_samples_per_treatment(treatment="control")

GROUPS = {
    "group1" : ["ChIP1", "ChIP2", "ChIP3"],
    "group2" : ["ChIP4", "ChIP5", "ChIP6"]
}                                           #I used this dictionnary to define the group of sample used in the multiBamSummary, might be improved a lot

##############
# Wildcards
##############
wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "L[0-9]+"

##############
# Desired output
##############

FASTQC_REPORTS  =     expand(RESULT_DIR + "fastqc/{sample}.html", sample=SAMPLES)
BAM_INDEX       =     expand(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam.bai", sample=SAMPLES)
BAM_RMDUP       =     expand(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam", sample=SAMPLES)
BIGWIG          =     expand(RESULT_DIR + "bigwig/{sample}.bw", sample=SAMPLES)
BAM_COMPARE     =     expand(RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw", zip, treatment = CASES, control = CONTROLS) #add zip function in the expand to compare respective treatment and control
BED_NARROW      =     expand(RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.narrowPeak", zip, treatment = CASES, control = CONTROLS)
BED_BROAD       =     expand(RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.broadPeak", zip, treatment = CASES, control = CONTROLS)
MULTIBAMSUMMARY =     RESULT_DIR + "multiBamSummary/MATRIX.npz"
PLOTCORRELATION =     RESULT_DIR + "plotCorrelation/MATRIX.png"
COMPUTEMATRIX   =     expand(RESULT_DIR + "computematrix/{treatment}_{control}.TSS.gz", treatment = CASES, control = CONTROLS)
HEATMAP         =     expand(RESULT_DIR + "heatmap/{treatment}_{control}.pdf", treatment = CASES, control = CONTROLS)
PLOTFINGERPRINT =     expand(RESULT_DIR + "plotFingerprint/{treatment}_vs_{control}.pdf", zip, treatment = CASES, control = CONTROLS)
PLOTPROFILE_PDF =     expand(RESULT_DIR + "plotProfile/{treatment}_{control}.pdf", treatment = CASES, control = CONTROLS)
PLOTPROFILE_BED =     expand(RESULT_DIR + "plotProfile/{treatment}_{control}.bed", treatment = CASES, control = CONTROLS)
MULTIQC         =     RESULT_DIR + "multiqc_report.html"

###############
# Final output
################
rule all:
    input:
        BAM_INDEX,
        BAM_RMDUP,
        FASTQC_REPORTS,
        #BEDGRAPH,
        BIGWIG,
        BED_NARROW,
        #BED_BROAD
        MULTIBAMSUMMARY,
        PLOTCORRELATION,
        COMPUTEMATRIX,
        HEATMAP,
        PLOTFINGERPRINT,
        PLOTPROFILE_PDF,
        PLOTPROFILE_BED,
        #MULTIQC
    message: "ChIP-seq pipeline succesfully run."		#finger crossed to see this message!

    shell:
      "#rm -rf {WORKING_DIR}"

###############
# Rules
###############

include : "rules/external_data.smk"
include : 'rules/pre_processing.smk'
include : "rules/macs2_peak_calling.smk"
include : "rules/deeptools_post_processing.smk"

rule multiqc:
    input:
        expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"}),
        expand(RESULT_DIR + "bed/{sample}_peaks.xls", sample= SAMPLES)
    output:
        "qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.27.1/bio/multiqc"
