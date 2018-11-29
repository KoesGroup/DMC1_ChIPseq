rule trimmomatic_se:
    input:
        get_fastq
    output:
        WORKING_DIR + "trimmed/{sample}.fastq.gz"
    message: "Trimming single-end {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic_se/{sample}.log"
    params :
        trimmer = ["TRAILING:3"],
        extra = "",
        # seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        # palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        # simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        # LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        # TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        # windowSize =                str(config['trimmomatic']['windowSize']),
        # avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        # minReadLen =                str(config['trimmomatic']['minReadLength']),
        # phred = 		            str(config["trimmomatic"]["phred"]),
        adapters =                  config["adapters"]
    threads: 10
    # conda:
    #     "../envs/trimmomatic_env.yaml"
    # shell:
    #     "trimmomatic SE {params.phred} -threads {threads} "
    #     "{input} "
    #     "{output} "
    #     "ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
    #     "LEADING:{params.LeadMinTrimQual} "
    #     "TRAILING:{params.TrailMinTrimQual} "
    #     "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
    #     "MINLEN:{params.minReadLen}"
    wrapper:
        "0.27.1/bio/trimmomatic/se"

rule trimmomatic_pe:
    input:
        get_fastq,

    output:
        forward_reads  = WORKING_DIR + "trimmed/{sample}.1.fastq.gz",
        reverse_reads  = WORKING_DIR + "trimmed/{sample}.2.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}.1.unpaired.fastq.gz"),
        reverseUnpaired  = temp(WORKING_DIR + "trimmed/{sample}.2.unpaired.fastq.gz")
    message: "Trimming paired-end {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic_pe/{sample}.log"
    params :
        trimmer = ["TRAILING:3"],
        extra = "",
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"]),
        adapters =                  config["adapters"]
    threads: 10
    conda:
        "../envs/trimmomatic.yaml"
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{params.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} &>{log}"
    # wrapper:
    #     "0.27.1/bio/trimmomatic/pe"

rule fastqc:
    input:
        WORKING_DIR + "trimmed/{sample}.fastq.gz",
    output:
        html = RESULT_DIR + "fastqc/{sample}.fastqc.html",
        rev = RESULT_DIR + "fastqc/{sample}.fastqc.zip"
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastqc.log"
    params:
        ""
    message:
        "Quality check of trimmed {wildcards.sample} sample with FASTQC"
    wrapper:
        "0.27.1/bio/fastqc"

rule index:
    input:
        WORKING_DIR + "genome.fasta"
    output:
        [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "genome.rev.1.bt2",
        WORKING_DIR + "genome.rev.2.bt2"
    message:"Indexing Reference genome"
    params:
        WORKING_DIR + "genome"
    threads: 10
    conda:
        "../envs/samtools_bowtie.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        sample= get_trimmed_reads,
        # index           = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        WORKING_DIR + "mapped/{sample}.bam"
    message: "Mapping files {wildcards.sample} to Reference genome"
    params:
        # bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "genome",
        extra =""
    threads: 10
    log:
        RESULT_DIR + "logs/bowtie/{sample}.log"
    conda:
        "../envs/samtools_bowtie.yaml"
    # shell:
    #     "bowtie2 {params.bowtie} "
    #     "--threads {threads} "
    #     "-x {params.index} "
    #     "-1 {input.forward} -2 {input.reverse} "
    #     "-U {input.forwardUnpaired},{input.reverseUnpaired} "   # also takes the reads unpaired due to trimming
    #     "| samtools view -Sb - > {output} 2>{log}"                       # to get the output as a BAM file directly
    wrapper:
        "0.27.1/bio/bowtie2/align"

rule sort:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        temp(RESULT_DIR + "mapped/{sample}.sorted.bam")
    message:"Sorting {wildcards.sample} bam file"
    threads: 10
    log:
        RESULT_DIR + "logs/samtools/{sample}.sort.log"
    conda:
        "../envs/samtools.yaml"
    shell:"samtools sort -@ {threads} -o {output} {input} &>{log}"

rule rmdup:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.bam"
    output:
        bam = RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam",
        bai = RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam.bai"        #bai files required for the bigwig and bamCompare rules
    message: "Removing duplicate from file {wildcards.sample}"
    log:
        RESULT_DIR + "logs/samtools/{sample}.sorted.rmdup.bam.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools rmdup {input} {output.bam} &>{log}
        samtools index {output.bam}
        """
        #samtools manual says "This command is obsolete. Use markdup instead

# rule bedgraph:
#     input:
#         RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
#     output:
#         RESULT_DIR + "bedgraph/{sample}.sorted.rmdup.bedgraph"
#     params:
#         genome = WORKING_DIR + "genome"
#     message:
#         "Creation of {wildcards.sample} bedgraph file"
#     log:
#         RESULT_DIR + "logs/deeptools/{sample}.sorted.rmdup.bedgraph.log"
#     conda:
#         "../envs/bedtools_env.yaml"
#     shell:
#         "bedtools genomecov -bg -ibam {input} -g {params.genome} > {output}"
