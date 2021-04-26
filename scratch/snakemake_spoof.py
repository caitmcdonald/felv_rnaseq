# snakemake spoof
# do this in your python interpreter...

### create a generic object class (called Foo here) ###
class Foo(object):
    pass

# make a variable, wildcards, which is an object of that class
wildcards = Foo()
wildcards.sample_id = 'S01'

### test function definitions ###
from snakemake.io import load_configfile
config = load_configfile("config/config.yaml")

import pandas as pd

## path to sequencing data:
metadata = r"{run_dir}/{sample_file}".format(run_dir= config["run_dir"], sample_file= config["samples"]
)

## read in samples
samples = pd.read_table(metadata).set_index("sample_id", drop=False)

## test definitions
def get_raw_fq1(wildcards):
    """Get path to R1 fastq files"""
    return r"raw/{fq}".format(fq=samples.loc[wildcards.sample_id, "fastq1"])

get_raw_fq1(wildcards)

def get_fq_all(wildcards):
    import glob, sys
    trimmed = glob.glob('data/trimmed/*val*')

    filename = []
    for f in trimmed:
        no_ext = f.split('.')[0]
        filename.append(no_ext)

    return filename
