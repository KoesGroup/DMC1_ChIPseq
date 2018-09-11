# Snakemake file for ChIP-Seq PE analysis

###############
# Libraries
###############

import os
import pandas as pd
from snakemake.utils import validate, min_version
#############################################
# Configuration and sample sheets
#############################################

configfile: "configs/config_tomato_sub.yaml"

WORKING_DIR = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR = config["result_dir"]      # what you want to keep

GENOME_FASTA_URL = config["refs"]["genome_url"]
GENOME_FASTA_FILE = os.path.basename(config["refs"]["genome_url"])
TOTALCORES = 16                             #check this via 'grep -c processor /proc/cpuinfo'

###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_samples_per_treatment(input_df="units.tsv",colsamples="sample",coltreatment="condition",treatment="control"):
    """This function returns a list of samples that correspond to the same experimental condition"""
    df = pd.read_table(input_df)
    df = df.loc[df[coltreatment] == treatment]
    filtered_samples = df[colsamples].tolist()
    return filtered_samples

##############
# Samples and conditions
##############

#samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
#SAMPLES = list(set(samples.index.values))

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
SAMPLES = units.index.get_level_values('sample').unique().tolist()
UNITS = units.index.get_level_values('unit').unique().tolist()

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

FASTQC_REPORTS  =     expand(RESULT_DIR + "fastqc/{sample}_{unit}_{pair}_fastqc.zip", sample=SAMPLES,unit=UNITS, pair={"forward", "reverse"})
BAM_INDEX       =     expand(RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam.bai", sample=SAMPLES,unit=UNITS)
BAM_RMDUP       =     expand(RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam", sample=SAMPLES,unit=UNITS)
BEDGRAPH        =     expand(RESULT_DIR + "bedgraph/{sample}_{unit}.sorted.rmdup.bedgraph", sample=SAMPLES,unit=UNITS)
BIGWIG          =     expand(RESULT_DIR + "bigwig/{sample}_{unit}.bw", sample=SAMPLES,unit=UNITS)
BAM_COMPARE     =     expand(RESULT_DIR + "bamcompare/log2_{treatment}_{control}_{unit}.bamcompare.bw", zip, treatment = CASES, control = CONTROLS,unit=UNITS) #add zip function in the expand to compare respective treatment and control
BED_NARROW      =     expand(RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.narrowPeak", treatment = CASES, control = CONTROLS,unit=UNITS)
BED_BROAD       =     expand(RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.broadPeak", zip, treatment = CASES, control = CONTROLS,unit=UNITS)

###############
# Final output
################
rule all:
    input:
        BAM_INDEX,
        BAM_RMDUP,
        FASTQC_REPORTS,
        BEDGRAPH,
        BIGWIG,
        BAM_COMPARE,
        BED_NARROW,
        #BED_BROAD
    message: "ChIP-seq pipeline succesfully run."		#finger crossed to see this message!

    shell:"#rm -rf {WORKING_DIR}"

###############
# Rules
###############
rule get_genome_fasta:
    output:
        WORKING_DIR + "genome.fasta"
    message:"downloading {GENOME_FASTA_FILE} genomic fasta file"
    shell: "wget -O {output} {GENOME_FASTA_URL}"

rule trimmomatic:
    input:
        reads = get_fastq,
        adapters = config["adapters"]
    output:
        forward_reads = WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq.gz",
        reverse_reads = WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}_{unit}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{sample}_{unit}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}_{unit}.log"
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
        "envs/trimmomatic_env.yaml"
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
        fwd=WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq.gz",
        rev=WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq.gz"
    output:
        fwd=RESULT_DIR + "fastqc/{sample}_{unit}_forward_fastqc.zip",
        rev=RESULT_DIR + "fastqc/{sample}_{unit}_reverse_fastqc.zip"
    log:
        RESULT_DIR + "logs/fastqc/{sample}_{unit}.fastqc.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "---Quality check of trimmed {wildcards.sample} sample with FASTQC"
    conda:
        "envs/fastqc_env.yaml"
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} 2>{log}"

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
        "envs/samtools_bowtie_env.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        forward = WORKING_DIR + "trimmed/{sample}_{unit}_forward.fastq.gz",
        reverse = WORKING_DIR + "trimmed/{sample}_{unit}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_{unit}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_{unit}_reverse_unpaired.fastq.gz",
        index = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        temp(WORKING_DIR + "mapped/{sample}_{unit}.bam")
    message: "Mapping files"
    params:
        bowtie = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index = WORKING_DIR + "genome"
    threads: 10
    conda:
        "envs/samtools_bowtie_env.yaml"
    shell:
        "bowtie2 {params.bowtie} "
        "--threads {threads} "
        "-x {params.index} "
        "-1 {input.forward} -2 {input.reverse} "
        "-U {input.forwardUnpaired},{input.reverseUnpaired} "   # also takes the reads unpaired due to trimming
        "| samtools view -Sb - > {output}"                       # to get the output as a BAM file directly

rule sort:
    input:
        WORKING_DIR + "mapped/{sample}_{unit}.bam"
    output:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.bam"
    message:"sorting {wildcards.sample} bam file"
    threads: 10
    conda:
        "envs/samtools.yaml"
    shell:"samtools sort -@ {threads} -o {output} {input}"

rule rmdup:
    input:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.bam"
    output:
        bam = RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam",
        bai = RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam.bai"        #bai files required for the bigwig and bamCompare rules
    message: "Removing duplicate from file {wildcards.sample}"
    log:
        RESULT_DIR + "logs/samtools/{sample}_{unit}.sorted.rmdup.bam.log"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools rmdup {input} {output.bam}
        samtools index {output.bam}
        """
        #samtools manual says "This command is obsolete. Use markdup instead

rule bedgraph:
    input:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bedgraph/{sample}_{unit}.sorted.rmdup.bedgraph"
    params:
        genome = WORKING_DIR + "genome"
    message:
        "Creation of {wildcards.sample} bedgraph file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}_{unit}.sorted.rmdup.bedgraph.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bedtools genomecov -bg -ibam {input} -g {params.genome} > {output}"

rule bigwig:
    input:
        RESULT_DIR + "mapped/{sample}_{unit}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bigwig/{sample}_{unit}.bw"
    message:
        "Converting {wildcards.sample} bam into bigwig file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}_{unit}_bamtobigwig.log"
    params:
        EFFECTIVEGENOMESIZE = str(config["bamCoverage"]["params"]["EFFECTIVEGENOMESIZE"]) #take argument separated as a list separated with a space
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input} -o {output} --effectiveGenomeSize {params.EFFECTIVEGENOMESIZE} 2>{log}"

rule bamcompare:
    input:
        treatment = RESULT_DIR + "mapped/{treatment}_{unit}.sorted.rmdup.bam",              #input requires an indexed bam file
        control = RESULT_DIR + "mapped/{control}_{unit}.sorted.rmdup.bam"                   #input requires an indexed bam file
    output:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}_{unit}.bamcompare.bw"
    message:
        "Running bamCompare"
    log:
        RESULT_DIR + "logs/deeptools/log2_{treatment}_{control}_{unit}.bamcompare.bw.log"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCompare -b1 {input.treatment} -b2 {input.control} -o {output.bigwig}"

rule call_narrow_peaks:
    input:
        treatment = RESULT_DIR + "mapped/{treatment}_{unit}.sorted.rmdup.bam",
        control = RESULT_DIR + "mapped/{control}_{unit}.sorted.rmdup.bam"
    output:
        bed = RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.narrowPeak"
    message:
        "Calling narrowPeak"
    params:
        name = "{treatment}_vs_{control}_{unit}",        #this option will give the output name, has to be similar to the output
        format = str(config['macs2']['format']),
        genomesize = str(config['macs2']['genomesize']),
        qvalue = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_{unit}_peaks.narrowPeak.log"
    conda:
        "envs/macs2_env.yaml"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/
        """

rule call_broad_peaks:
    input:
        treatment = RESULT_DIR + "mapped/{treatment}_{unit}.sorted.rmdup.bam",
        control = RESULT_DIR + "mapped/{control}_{unit}.sorted.rmdup.bam"
    output:
        bed = RESULT_DIR + "bed/{treatment}_vs_{control}_{unit}_peaks.broadPeak"
    message:
        "Calling broadPeak"
    params:
        name = "{treatment}_vs_{control}_{unit}",
        format = str(config['macs2']['format']),
        genomesize = str(config['macs2']['genomesize']),
        qvalue = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_{unit}_peaks.broadPeak.log"
    conda:
        "envs/macs2_env.yaml"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} --broad --broad-cutoff 0.1 {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/
        """
