#!/bin/bash

#SBATCH --mem 92GB
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coby.mcdonald@colostate.edu
#SBATCH --output=logs/trimgalore/output-%j
#SBATCH --error=logs/trimgalore/error-%j

source ~/.bashrc
conda activate trim_galore

export WORK=/scratch/summit/camcd@colostate.edu/felv_rnaseq/data/raw
cd $WORK

for i in $(ls *R1_001.fastq.gz | sed 's/\R1_001.fastq.gz//'); do
    trim_galore --cores 20 --paired --retain_unpaired --phred33 --length 36 -q 5 --stringency 1 -e 0.1 -o ../../results/trimgalore ./$i\R1_001.fastq.gz ./$i\R2_001.fastq.gz;
done

# runtime: 2.5 hrs

# okay so what this script is doing: for a directory that has PE fastq.gz files named <whatever-the-name-is.R1_001.fastq.gz> and <whatever-the-name-is.R2_001.fastq.gz>, loop over each of those filenames and run trim_galore on each sample pair. the stuff at the top is all if you want to run it using SLURM, so if you don't you can just run it straight as a bash one-liner for-loop with nohup from whatever your directory of interest is. i don't remember, but i think you may need to mkdir the <trimmed> output directory before running trimgalore. i used all of the --length, -q, --stringency, and -e settings at default values.
