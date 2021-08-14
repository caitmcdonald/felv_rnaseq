#!/bin/bash

#SBATCH --mem 92GB
#SBATCH --time=23:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coby.mcdonald@colostate.edu
#SBATCH --output=output-%j
#SBATCH --error=error-%j


## Run OrthoFinder
source ~/.bashrc
conda activate orthofinder

export WORK=/scratch/summit/camcd@colostate.edu/felv_rnaseq/results/orthofinder
cd $WORK
# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/327/715/GCF_003327715.1_PumCon1.0/GCF_003327715.1_PumCon1.0_protein.faa.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_000181335.3_Felis_catus_9.0/GCF_000181335.3_Felis_catus_9.0_protein.faa.gz
# gunzip *.gz

## ID primary transcripts
# for f in *.faa ; do python primary_transcript.py $f ; done

## Run orthofinder on primary transcripts using Diamond database
orthofinder -S diamond -f primary_transcripts/ #default 16 threads
# orthofinder -t 20 -f primary_transcripts/
