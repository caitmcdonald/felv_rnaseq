#!/bin/bash

#SBATCH --mem 5GB
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coby.mcdonald@colostate.edu
#SBATCH --output=logs/qc/raw/output-%j
#SBATCH --error=logs/qc/raw/error-%j


## load envs
source ~/.bashrc
conda activate qctrim

## run fastqc
fastqc -t 19 data/raw/*.fastq.gz -o results/fastqc/raw

## run multiqc
multiqc results/fastqc/raw -o results/multiqc/raw
