#!/bin/bash

#SBATCH --mem 74GB
#SBATCH --time=23:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cait.mcdonald@colostate.edu
#SBATCH --output=output-%j
#SBATCH --error=error-%j


## Run OrthoFinder
source ~/.bashrc
conda activate orthofinder

cd /scratch/summit/camcd@colostate.edu/felv_rnaseq/resources/orthofinder_tutorial/

# wget http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz

# wget http://ftp.ensembl.org/pub/current_fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz

# for f in *fa ; do python primary_transcript.py $f ; done

orthofinder -S diamond -f primary_transcripts/ #default 16 threads
