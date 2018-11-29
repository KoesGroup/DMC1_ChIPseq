# Snakemake file for ChIP-Seq PE analysis

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

WORKING_DIR         = config["temp_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR          = config["result_dir"]      # what you want to keep

GENOME_FASTA_URL    = config["refs"]["genome_url"]
GENOME_FASTA_FILE   = os.path.basename(config["refs"]["genome_url"])
TOTALCORES          = config["cores"] 

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

##############
# Samples and conditions
##############

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

CASES = get_samples_per_treatment(treatment="treatment")
CONTROLS = get_samples_per_treatment(treatment="control")

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

###############
# Final output
################
rule all:
    input:
        FASTQC_REPORTS,
        BIGWIG,
        BED_NARROW,
        BED_BROAD
    message: "ChIP-seq pipeline succesfully run."		#finger crossed to see this message!

    shell:"rm -rf {WORKING_DIR}"

###############
# Rules
###############
rule get_genome_fasta:
    output:
        WORKING_DIR + "genome.fasta"
    message:"downloading {GENOME_FASTA_FILE} genomic fasta file"
    conda: 
        "envs/wget.yaml"
    shell: "wget -O {output} {GENOME_FASTA_URL}"

rule trimmomatic:
    input:
        reads = get_fastq,
        adapters = config["trimmomatic"]["adapters"]
    output:
        forward_reads   = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse_reads   = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}.log"
    params :
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.reads} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} 2>{log}"

rule fastqc:
    input:
        fwd = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        rev = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz"
    output:
        RESULT_DIR + "fastqc/{sample}.html"
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastp.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "Quality check of trimmed {wildcards.sample} sample with fastp"
    conda:
        "envs/fastp.yaml"
    shell:
        "fastp -i {input.fwd} -I {input.rev} -h {output} 2>{log}"

rule index:
    input:
        WORKING_DIR + "genome.fasta"
    output:
        [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "genome.rev.1.bt2",
        WORKING_DIR + "genome.rev.2.bt2"
    message:"indexing genome"
    params:
        WORKING_DIR + "genome"
    threads: 10
    conda:
        "envs/bowtie2.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        forward         = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse         = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz",
        index           = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        temp(WORKING_DIR + "mapped/{sample}.bam")
    message: "Mapping files"
    params:
        bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "genome"
    threads: 10
    conda:
        "envs/samtools_bowtie.yaml"
    shell:
        "bowtie2 {params.bowtie} "
        "--threads {threads} "
        "-x {params.index} "
        "-1 {input.forward} -2 {input.reverse} "
        "-U {input.forwardUnpaired},{input.reverseUnpaired} "    # also takes the reads unpaired due to trimming
        "| samtools view -Sb - > {output}"                       # to get the output as a BAM file directly

rule sort:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        temp(RESULT_DIR + "mapped/{sample}.sorted.bam")
    message:"sorting {wildcards.sample} bam file"
    threads: 10
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule rmdup:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.bam"
    output:
        bam = temp(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"),
        bai = temp(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam.bai")        #bai files required for the bigwig and bamCompare rules
    message: "Removing duplicate from file {wildcards.sample} using samtools rmdup"
    log:
        RESULT_DIR + "logs/samtools/{sample}.sorted.rmdup.bam.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools rmdup {input} {output.bam} 2>{log}
        samtools index {output.bam}
        """

rule bedgraph:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bedgraph/{sample}.sorted.rmdup.bedgraph"
    params:
        genome = WORKING_DIR + "genome"
    message:
        "Creation of {wildcards.sample} bedgraph file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}.sorted.rmdup.bedgraph.log"
    conda:
        "envs/bedtools.yaml"
    shell:
        "bedtools genomecov -bg -ibam {input} -g {params.genome} > {output}"

rule bigwig:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bigwig/{sample}.bw"
    message:
        "Converting {wildcards.sample} bam into bigwig file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}_bamtobigwig.log"
    params:
        EFFECTIVEGENOMESIZE = str(config["bamCoverage"]["params"]["EFFECTIVEGENOMESIZE"]), #take argument separated as a list separated with a space
        EXTENDREADS         = str(config["bamCoverage"]["params"]["EXTENDREADS"])
    conda:
        "envs/deeptools.yaml"
    shell:
        "samtools index {input};"
        "bamCoverage --bam {input} -o {output} --effectiveGenomeSize {params.EFFECTIVEGENOMESIZE} --extendReads {params.EXTENDREADS} 2>{log}"

rule bamcompare:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",              #input requires an indexed bam file
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"                   #input requires an indexed bam file
    output:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw"
    message:
        "Running bamCompare on {input.control} and {input.treatment} bam files"
    log:
        RESULT_DIR + "logs/deeptools/log2_{treatment}_{control}.bamcompare.bw.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCompare -b1 {input.treatment} -b2 {input.control} -o {output.bigwig} 2>{log}"

rule call_narrow_peaks:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.narrowPeak"
    message:
        "Calling narrow peaks with macs2"
    params:
        name        = "{treatment}_vs_{control}",        # this option will give the output name, has to be similar to the output
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_peaks.narrowPeak.log"
    conda:
        "envs/macs2.yaml"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/ 2>{log}
        """

rule call_broad_peaks:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.broadPeak"
    message:
        "Calling broad peaks with macs2"
    params:
        name        = "{treatment}_vs_{control}",
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue']),
        resultdir   = RESULT_DIR + "bed/"  
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_peaks.broadPeak.log"
    conda:
        "envs/macs2.yaml"
    shell:
        "macs2 callpeak -t {input.treatment} -c {input.control} " 
        "--params {params.format} "
        "--broad --broad-cutoff 0.1 "
        "--gsize {params.genomesize} "
        "--name {params.name} "
        "--nomodel "
        "--bdg "
        "-q {params.qvalue} "
        "--outdir {params.outdir} "
        "2>{log}"



