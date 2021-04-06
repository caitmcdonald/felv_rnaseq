### Approach 1: use Python to mess with file names to make them easier for wildcards:

# import glob, sys #following Eric A--omg PYTHON exciting!!!
#
# fullnames= glob.glob('/data/*.fastq.gz')
#
# # get sample name without .fastq.gz
# sample= []
# for name in fullnames:
#     base= name.split('.')[0]
#     sample.append(base)
#
# # get base name without R1 or R2
# basename= []
# for read in sample:
#     noR= read.split('_')[0] + '_' + read.split('_')[1]
#     basename.append(noR)

### Approach 2: import metadata as pandas df
import pandas as pd
#
metadata = pd.read_table("felv_metadata_long.txt").set_index("sample_id", drop=False)
# fastq_list = metadata['fastq'].tolist() #stores column as list of attributes
# Integrating pandas df with snakemake supposedly:
# input:
    # lambda wildcards, output: metadata.fastq[wildcards.sample_id]
class Foo(object):
    pass

wildcards= Foo()

wildcards.sample_id= metadata['sample_id'].tolist()

### Approach 3: just list all the samples in a config file

configfile: "config.yaml"

# testing:
rule fastqc_raw:
    input:
        # "data/samples/{sample}.fastq.gz"
        # get_fastqc_input_fastqs
        # lambda wildcards: config["samples"][wildcards.sample]
        lambda wildcards: metadata.fastq[wildcards.sample_id]
    output:
        "results/fastqc/raw/{sample_id}.html",
        "results/fastqc/raw/{sample_id}.zip"
    conda:
        "qctrim.yaml"
    log:
        "logs/fastqc/raw/{sample_id}.log"
    threads: 1
    shell:
        "(fastqc {input}) 2> {log}"


# don't totally understand this rule yet...cannot have wildcards within target rule
# rule all:
#     input:
#         #the last output
#         "results/multiqc/raw/multiqc_report.html",
#         directory("results/multiqc/raw/multiqc_data")
#
# rule fastqc_raw:
#     input:
#         # "/data/{sample}.fastq.gz"
#         # get_fastqc_input_fastqs
#         lambda wildcards: config["samples"][wildcards.sample]
#     output:
#         "results/fastqc/raw/{sample}.html",
#         "results/fastqc/raw/{sample}.zip"
#     conda:
#         "qctrim.yaml"
#     log:
#         "logs/fastqc/raw/{sample}.log"
#     threads: 1
#     shell:
#         "(fastqc {input}) 2> {log}"
#
# rule multiqc_raw:
#     input:
#         "results/fastqc/raw"
#     output:
#         "results/multiqc/raw/multiqc_report.html",
#         directory("results/multiqc/raw/multiqc_data")
#     shell:
#         "multiqc {input}"

# rule cutadapt:
#     input:
#         "/data/{sample}.fastq.gz"
#     output:
#     shell:
#
# rule fastqc_trim:
#     input:
#         "/results/trimming/{basename}_val_1.fastq.gz",
#         "/results/trimming/{basename}_val_2.fastq.gz"
#     output:
#         "/results/fastqc/trimmed/{sample}.trimmed_fastqc.zip",
#         "/results/fastqc/trimmed/{sample}.trimmed_fastqc.html"
#     shell:
#         "fastqc {input}"
#
# rule multiqc_trim:
#     input:
#         "/results/fastqc/trimmed/"
#     output:
#         "/results/multiqc/trimmed/multiqc_report.html",
#         #directory("/results/multiqc/trimmed/multiqc_data")
#     shell:
#         "multiqc {input}"
