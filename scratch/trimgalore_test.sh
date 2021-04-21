#!/bin/bash

#SBATCH --mem 92GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cait.mcdonald@colostate.edu
#SBATCH --output=output-%j
#SBATCH --error=error-%j

source ~/.bashrc
conda activate trim_galore

trim_galore --cores 20 --paired --retain_unpaired --phred33 --length 36 -q 5 --stringency 1 -e 0.1 -o data/trimmed data/raw/4438_S1_L002_R1_001.fastq.gz data/raw/4438_S1_L002_R2_001.fastq.gz
