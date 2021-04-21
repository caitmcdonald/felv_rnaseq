# snakemake spoof
# do this in your python interpreter...

### create a generic object class (called Foo here) ###
class Foo(object):
    pass

# make a variable, wildcards, which is an object of that class
wildcards = Foo()
wildcards.library = '3480_ATCACG_S83_L008'
wildcards.library


### test function definitions ###
from snakemake.io import load_configfile
config = load_configfile("config/config_test.yaml")

import pandas as pd

## path to sequencing data:
seq_dat = r"{run_dir}/{sample_file}".format(run_dir= config["run_dir"], sample_file= config["samples"]
)
barcode_dat = r"{run_dir}/{barcode_file}".format(run_dir= config["run_dir"], barcode_file= config["barcodes"]
)

## read in samples
samples = pd.read_table(seq_dat).set_index("library", drop=False)
animals = pd.read_table(barcode_dat)
#
animal_id = (animals.groupby(['library', 'path_fq1', 'path_fq2', 'barcode_file', 'lib_dir'])
        .animal_id.apply(lambda x: x.tolist())
        .reset_index())
barcode = (animals.groupby(['library'])
        .barcode.apply(lambda x: x.tolist())
        .reset_index())
collection_date = (animals.groupby(['library'])
        .collection_date.apply(lambda x: x.tolist())
        .reset_index())
animal_list = animal_id.merge(barcode, how='left')
animal_list = animal_list.merge(collection_date, how='left')
animals = animal_list.set_index("library", drop=False)

## test definitions

### I need to figure out functions to take the library id (e.g. sample2_S2_L005) and return a) the correct barcode file...should I do this as a dictionary???
def get_barcodes(wildcards):
    """Get path to barcodes .txt"""
    return r"data/raw/{lib_dir}/{barcode_file}".format(lib_dir=animals.loc[wildcards.library, "lib_dir"], barcode_file=animals.loc[wildcards.library, "barcode_file"])
    # return r"data/raw/{lib_dir}/{barcode_file}".format(lib_dir=samples.loc[wildcards.library, "lib_dir"], barcode_file=samples.loc[wildcards.library, "barcode_file"])
get_barcodes(wildcards)

def get_animals(wildcards):
    """Get all animals multiplexed per library"""
    return animals.loc[wildcards.library, "animal_id"]
get_animals(wildcards)
