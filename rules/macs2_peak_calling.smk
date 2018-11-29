rule call_narrow_peaks:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.narrowPeak"
    message:
        "Calling narrowPeak for {wildcards.treatment} and {wildcards.control}"
    params:
        name        = "{treatment}_vs_{control}",        #this option will give the output name, has to be similar to the output
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_peaks.narrowPeak.log"
    conda:
        "../envs/macs2.yaml"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/ &>{log}
        """

rule call_broad_peaks:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.sorted.rmdup.bam",
        control     = RESULT_DIR + "mapped/{control}.sorted.rmdup.bam"
    output:
        RESULT_DIR + "bed/{treatment}_vs_{control}_peaks.broadPeak"
    message:
        "Calling broadPeak for {wildcards.treatment} and {wildcards.control}"
    params:
        name        = "{treatment}_vs_{control}",
        format      = str(config['macs2']['format']),
        genomesize  = str(config['macs2']['genomesize']),
        qvalue      = str(config['macs2']['qvalue'])
    log:
        RESULT_DIR + "logs/macs2/{treatment}_vs_{control}_peaks.broadPeak.log"
    conda:
        "../envs/macs2.yaml"
    shell:
        """
        macs2 callpeak -t {input.treatment} -c {input.control} {params.format} --broad --broad-cutoff 0.1 {params.genomesize} --name {params.name} --nomodel --bdg -q {params.qvalue} --outdir results/bed/ &>{log}
        """
