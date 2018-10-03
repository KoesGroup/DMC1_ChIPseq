# ChIP_seq_Snakemake
A snakemake pipeline for the analysis of ChIP-seq data

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Miniconda](https://img.shields.io/badge/miniconda-blue.svg)](https://conda.io/miniconda)

# Aim
Snakemake pipeline made for reproducible analysis of paired-end Illumina ChIP-seq data. The desired output of this pipeline are:
- fastqc zip and html files
- bigWig files (including bamCompare rule)
- bed files from narrow and broad peaks
- bigWig files 

# Content of the repository

- **Snakefile** containing the targeted output and the rules to generate them from the input files.

- **config/** , folder containing the configuration files making the Snakefile adaptable to any input files, genome and parameter for the rules. Adapt the config file and its reference in the Snakefile.

- **Fastq/**, folder containing subsetted paired-end fastq files used to test locally the pipeline. Generated using [Seqtk](https://github.com/lh3/seqtk): `seqtk sample -s100 read1.fq 5000 > sub1.fqseqtk sample -s100 read2.fq 5000 > sub2.fq`. RAW fastq or fastq.gz files should be placed here before running the pipeline.

- **envs/**, folder containing the environment needed for the Snakefile to run. To use Snakemake, it is required to create and activate an environment containing snakemake (here : envs/global_env.yaml )

- **units.tsv**, is a tab separated value files containing information about the experiment name, the condition of the experiment (control or treatment) and the path to the fastq files relative to the **Snakefile**. **Change this file according to your samples.**


# Usage

## Conda environment

**First**, you need to create an environment for the use of `Snakemake` with [Conda package manager](https://conda.io/docs/using/envs.html).
1. Create a virtual environment named "chipseq" from the `global_env.yaml` file with the following command: `conda env create --name chipseq --file ~/envs/global_env.yaml`
2. Then, activate this virtual environment with `source activate chipseq`

The `Snakefile` will then take care of installing and loading the packages and softwares required by each step of the pipeline.

## Configuration file
The `config.yaml` file specifies the sample list, the genomic reference fasta file to use, the directories to use, etc. This file is then used to build parameters in the main `Snakefile`.

## Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order.
It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/).

## Samples
Samples are listed in the `units.tsv` file and will be used by the Snakefile automatically. Change in the 'units.tsv', the columns `fq1`, `fq2` and `conditions` according to your needs. `fq1` and `fq2` requires the path to the samples names relative to the Snakefile.

## Dry run
Use the command `snakemake -np` to perform a dry run that prints out the rules and commands.

## Real run
Simply type `snakemake --use-conda` and provide the number of cores with `--cores 10` for ten cores for instance.
For cluster execution, please refer to the [Snakemake reference](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution).
**Please pay attention** to add  `--use-conda`, it is required for the installation and loading of the dependencies used by the rules of the pipeline.
To run the pipeline, from the folder containing the Snakefile run the command : `snakemake --use-conda --core 10`.

# Main outputs
The main output are in this release (1.0) :

- **fastqc** : provides quality check analysis for `fq1` and `fq2` of each ChIP-seq experiment, the outputs files are `html` which can be open and viewed in any web browser. It gives a quick impression of wheter your data have any problem which you should be aware before doing any further analysis.

- **bed** : MACS2 peak calling function generates a list of files whose details can be found [here](https://github.com/taoliu/MACS). The **NAME_peaks.narrowPeak** can be used to visualize the location of peaks on the IGV.  

- **bigwig files** : Files can be used to display the read coverage over the genome on the IGV. They are faster to load on IGV than the bedgraph files as well. To combine with the **NAME_peaks.narrowPeak** to have a estimation of the read coverage at the peak location.

Optionals outputs of the pipelines are **bamCompare**, **bedgraph** and **bed files for broad peaks calling**, but they require to change the `Snakefile`.

## DAG of the pipeline
![](dag.png)
