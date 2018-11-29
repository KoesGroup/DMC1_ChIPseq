rule get_genome_fasta:
    output:
        WORKING_DIR + "genome.fasta"
    message:"Downloading {GENOME_FASTA_FILE} genomic fasta file"
    shell: "wget -O {output} {GENOME_FASTA_URL}"


rule get_gff:
    output:
        WORKING_DIR + "gene_model.gff"
    message:"Downloading {GFF_FILE} genomic fasta file"
    shell: "wget -O {output} {GFF_URL}"

rule gff_to_gtf:
    input:
        WORKING_DIR + "gene_model.gff"
    message:"Converting GFF file to GFF"    
    output:
        WORKING_DIR + "gene_model.gtf"
    shell:
        "python scripts/gff_to_gtf.py {input} {output}"    