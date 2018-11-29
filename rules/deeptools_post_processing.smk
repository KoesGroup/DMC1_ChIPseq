rule bamCoverage:
    input:
        RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bigwig/{sample}.bw"
    message:
        "Converting {wildcards.sample} bam into bigwig file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}_bamtobigwig.log"
    params:
        EFFECTIVEGENOMESIZE     = str(config["bamCoverage"]["params"]["EFFECTIVEGENOMESIZE"]), #take argument separated as a list separated with a space
        EXTENDREADS             = str(config["bamCoverage"]["params"]["EXTENDREADS"]),
        binSize                 = str(config['bamCoverage']["params"]['binSize']),
        normalizeUsing          = str(config['bamCoverage']["params"]['normalizeUsing']),
        ignoreForNormalization  = str(config['bamCoverage']["params"]['ignoreForNormalization']),
        smoothLength            = str(config['bamCoverage']["params"]['smoothLength'])
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input} \
        -o {output} \
        --effectiveGenomeSize {params.EFFECTIVEGENOMESIZE} \
        --extendReads {params.EXTENDREADS} \
        --binSize {params.binSize} \
        --smoothLength {params.smoothLength} \
        --ignoreForNormalization {params.ignoreForNormalization} \
        &>{log}"

rule bamcompare:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",              #input requires an indexed bam file
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"                   #input requires an indexed bam file
    output:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw"
    message:
        "Running bamCompare for {wildcards.treatment} and {wildcards.control}"
    log:
        RESULT_DIR + "logs/deeptools/log2_{treatment}_{control}.bamcompare.bw.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binSize             = str(config['bamcompare']['binSize']),
        normalizeUsing      = str(config['bamcompare']['normalizeUsing']),
        EFFECTIVEGENOMESIZE = str(config["bamcompare"]["EFFECTIVEGENOMESIZE"]),
        operation           = str(config['bamcompare']['operation']),
        smoothLength        = str(config['bamcompare']['smoothLength']),
        ignoreForNormalization = str(config['bamcompare']['ignoreForNormalization']),
        scaleFactorsMethod  = str(config['bamcompare']['scaleFactorsMethod'])
    shell:
        "bamCompare -b1 {input.treatment} \
        -b2 {input.control}  \
        --binSize {params.binSize} \
        -o {output.bigwig} \
        --normalizeUsing {params.normalizeUsing} \
        --operation {params.operation} \
        --smoothLength {params.smoothLength} \
        --ignoreForNormalization {params.ignoreForNormalization} \
        --scaleFactorsMethod {params.scaleFactorsMethod} \
        &>{log}"

rule multiBamSummary:
    input:
        lambda wildcards: expand(RESULT_DIR + "mapped/{sample}.sorted.rmdup.bam", sample = SAMPLES)
    output:
        RESULT_DIR + "multiBamSummary/MATRIX.npz"
    message:
        "Computing the read coverage into a numpy array "
    threads: 10
    params:
        binSize     = str(config['multiBamSummary']['binSize'])
    log:
        RESULT_DIR + "logs/deeptools/multibamsummary/MATRIX.log"
    shell:
        "multiBamSummary bins \
        --bamfiles {input} \
        --numberOfProcessors {threads}\
        --binSize {params.binSize} \
        --centerReads \
        --extendReads \
        -o {output} \
        2> {log}"


rule plotCorrelation:
    input:
        RESULT_DIR + "multiBamSummary/MATRIX.npz"
    output:
        RESULT_DIR + "plotCorrelation/MATRIX.png"
    log:
        RESULT_DIR + "logs/deeptools/plotcorrelation/MATRIX.log"
    params:
        corMethod  = str(config['plotCorrelation']['corMethod']),
        whatToPlot = str(config['plotCorrelation']['whatToPlot']),
        color      = str(config['plotCorrelation']['color'])
    conda:
        "../envs/deeptools.yaml"
    message:
        "Preparing the correlation plot between all samples"
    shell:
        "plotCorrelation \
                    --corData {input} \
                    --corMethod {params.corMethod} \
                    --whatToPlot {params.whatToPlot} \
                    --skipZeros \
                    --colorMap {params.color} \
                    --plotFile {output} \
                    --plotNumbers \
                    2> {log}"

#--plotTitle 'Pearson Correlation of {params.title} coverage' \


rule computeMatrix:
    input:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw",
        bed    = WORKING_DIR + "gene_model.gtf"
    output:
        RESULT_DIR + "computematrix/{treatment}_{control}.TSS.gz"
    threads: 10
    params:
        binSize = str(config['computeMatrix']['binSize']),
        upstream    = str(config['computeMatrix']['upstream']),
        downstream = str(config['computeMatrix']['downstream'])
    conda:
        "../envs/deeptools.yaml"
    log:
        RESULT_DIR + "logs/deeptools/computematrix/{treatment}_{control}.log"
    message:
        "Computing matrix for {input.bigwig} with {params.binSize} bp windows and {params.upstream} bp around TSS"        
    shell:
        "computeMatrix \
        reference-point \
        --referencePoint TSS \
        -S {input.bigwig} \
        -R {input.bed} \
        --afterRegionStartLength {params.upstream} \
        --beforeRegionStartLength {params.downstream} \
        --numberOfProcessors {threads} \
        --binSize {params.binSize} \
        -o {output} \
        2> {log}"

rule plotHeatmap:
    input:
        RESULT_DIR + "computematrix/{treatment}_{control}.TSS.gz"
    output:
        RESULT_DIR + "heatmap/{treatment}_{control}.pdf"
    params:
        kmeans = str(config['plotHeatmap']['kmeans']),
        color  = str(config['plotHeatmap']['color']),
        plot   = str(config['plotHeatmap']['plot']),
        cluster = "{treatment}_vs_{control}.bed"
    conda:
        "../envs/deeptools.yaml"
    log:
        RESULT_DIR + "logs/deeptools/plotHeatmap/{treatment}_{control}.log"
    message:
        "Preparing Heatmaps..."        
    shell:
        "plotHeatmap \
        --matrixFile {input} \
        --outFileName {output} \
        --kmeans {params.kmeans} \
        --colorMap {params.color} \
        --legendLocation best \
        --outFileSortedRegions {params.cluster}"

rule plotFingerprint:
    input:
        treatment = expand(RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam", treatment = CASES),
        control   = expand(RESULT_DIR + "mapped/{control}.sorted.rmdup.bam", control = CONTROLS)
    output:
        pdf = RESULT_DIR + "plotFingerprint/{treatment}_vs_{control}.pdf"
    params:
        EXTENDREADS  = str(config["bamCoverage"]["params"]["EXTENDREADS"]),
        binSize      = str(config['bamCoverage']["params"]['binSize'])
    conda:
        "../envs/deeptools.yaml"
    log:
        RESULT_DIR + "logs/deeptools/plotFingerprint/{treatment}_vs_{control}.log"
    message:
        "Preparing deeptools plotFingerprint"        
    shell:
        "plotFingerprint \
        -b {input.treatment} {input.control} \
        --extendReads {params.EXTENDREADS} \
        --binSize {params.binSize} \
        --plotFile {output}"

rule plotProfile:
    input:
        RESULT_DIR + "computematrix/{treatment}_{control}.TSS.gz"
    output:
        pdf = RESULT_DIR + "plotProfile/{treatment}_{control}.pdf",
        bed = RESULT_DIR + "plotProfile/{treatment}_{control}.bed"
    params:
        kmeans      = str(config['plotProfile']['kmeans']),
        startLabel  = str(config['plotProfile']['startLabel']),
        endLabel    = str(config['plotProfile']['endLabel'])
    conda:
        "../envs/deeptools.yaml"
    log:
        RESULT_DIR + "logs/deeptools/plotProfile/{treatment}_{control}.log"
    message:
        "Preparing deeptools plotProfile"        
    shell:
        "plotProfile \
        --matrixFile {input} \
        --outFileName {output.pdf} \
        --outFileSortedRegions {output.bed} \
        --kmeans {params.kmeans} \
        --startLabel {params.startLabel} \
        --endLabel {params.endLabel}"