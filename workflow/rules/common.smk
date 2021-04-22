import pandas as pd

## config file:
configfile: "config/config.yaml"

## path to sample metadata:
metadata = r"{run_dir}/{sample_file}".format(run_dir= config["run_dir"], sample_file= config["samples"]
)

## read in samples
samples = pd.read_table(metadata).set_index("sample_id", drop=False)

## functions for input wildcards

# get read1 fastq files

def get_raw_fq1(wildcards):
    """Get path to R1 fastq files"""
    return r"raw/{fq}".format(fq=samples.loc[wildcards.sample_id, "fastq1"])

# get read2 fastq files
def get_raw_fq2(wildcards):
    """Get path to R2 fastq files"""
    return r"raw/{fq}".format(fq=samples.loc[wildcards.sample_id, "fastq2"])

# get trimmed read1 fastq files
def get_trimmed_fq1(wildcards):
    """Get path to R1 fastq files"""
    return r"{run_dir}/trimmed/{fq}".format(run_dir=config["run_dir"], fq=samples.loc[wildcards.sample_id, "fastq1"])

# get trimmed read2 fastq files
def get_trimmed_fq2(wildcards):
    """Get path to R2 fastq files"""
    return r"{run_dir}/trimmed/{fq}".format(run_dir=config["run_dir"], fq=samples.loc[wildcards.sample_id, "fastq2"])
