#!/bin/bash

#SBATCH --mem 5GB
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coby.mcdonald@colostate.edu
#SBATCH --output=logs/qc/trimmed/output-%j
#SBATCH --error=logs/qc/trimmed/error-%j


## load envs
source ~/.bashrc
conda activate qctrim

## run fastqc
fastqc -t 12 results/trimgalore/*val_1.fq.gz -o results/fastqc/trimmed
fastqc -t 12 results/trimgalore/*val_2.fq.gz -o results/fastqc/trimmed

## run multiqc
multiqc results/fastqc/trimmed -o results/multiqc/trimmed
