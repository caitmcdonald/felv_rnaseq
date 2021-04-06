#!/bin/bash

#SBATCH --mem-per-cpu 800
#SBATCH --time=02:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cait.mcdonald@colostate.edu
#SBATCH --output=output-%j
#SBATCH --error=error-%j


## load envs
source ~/.bashrc
conda activate qctrim

## run fastqc
fastqc -t 19 data/samples/*.fastq.gz -o results/fastqc/raw

## run multiqc
multiqc results/fastqc/raw -o results/multiqc/raw -n results/multiqc/raw/multiqc_report.html
