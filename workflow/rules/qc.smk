import pandas as pd

## config:
configfile: config/config.yaml

## path to sample metadata:
metadata = r"{run_dir}/{sample_file}".format(run_dir= config["run_dir"], sample_file= config["samples"]
)

# met_long = pd.read_table("data/metadata/felv_metadata_long.txt").set_index("sample_id", drop=False)
samples = pd.read_table(metadata).set_index("sample_id", drop=False)

## now need to transform table to long format in pandas:


## functions for input wildcards

def get_raw_fastqs(wildcards):
    """Get path to raw fastq files"""
    return r"{run_dir}/raw/{fq}"

def get_input_fastqs(wildcards):
    return "data/samples/" + met_long.loc[wildcards.sample_id, "fastq"]

# don't totally understand this rule yet...cannot have wildcards within target rule
# rule all:
#     input:
#         #the last output
#         # directory("results/multiqc/raw/multiqc_data")
#         "results/multiqc/raw/multiqc_data"

## rules

rule fastqc_raw:
    input:
        get_input_fastqs
    output:
        "results/fastqc/raw/{sample_id}_fastqc.html",
        "results/fastqc/raw/{sample_id}_fastqc.zip"
    conda:
        "envs/fastqc.yaml" #note: don't list workflow parent dir as snakemake recognizes it automatically
    log:
        "logs/fastqc/raw/{sample_id}.log"
    threads: 19
    shell:
        "echo (fastqc -t {threads} {input}) 2> {log}"

rule multiqc_raw:
    input:
        "results/fastqc/raw"
    output:
        "results/multiqc/raw/multiqc_report.html",
        direct=directory("results/multiqc/raw")
    conda:
        "envs/multiqc.yaml"
    shell:
        "echo multiqc {input} -o {output.direct}"
