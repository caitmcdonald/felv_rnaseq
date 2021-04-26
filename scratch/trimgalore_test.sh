#!/bin/bash

#SBATCH --mem 92GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cait.mcdonald@colostate.edu
#SBATCH --output=output-%j
#SBATCH --error=error-%j

source ~/.bashrc
conda activate trim_galore

cd data/raw/

for i in $(ls *R1_001.fastq.gz | sed 's/\R1_001.fastq.gz//'); do
    trim_galore --cores 20 --paired --retain_unpaired --phred33 --length 36 -q 5 --stringency 1 -e 0.1 -o ../trimmed ./$i\R1_001.fastq.gz ./$i\R2_001.fastq.gz;
done
