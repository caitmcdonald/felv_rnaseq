#!/bin/bash

#SBATCH --mem 110GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cait.mcdonald@colostate.edu
#SBATCH --output=output-%j
#SBATCH --error=error-%j

## Get genome and gtf
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_000181335.3_Felis_catus_9.0/GCF_000181335.3_Felis_catus_9.0_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_000181335.3_Felis_catus_9.0/GCF_000181335.3_Felis_catus_9.0_genomic.gtf.gz
# gunzip *.gz


## Run STAR with quant mode
source ~/.bashrc
conda activate star_quant

export WORK=/scratch/summit/camcd@colostate.edu/felv_rnaseq
cd $WORK
mkdir resources/genome_star
export GENOMEDIR=resources/genome_star

## Build genome index
STAR --runThreadN 22 --runMode genomeGenerate --genomeDir $GENOMEDIR --genomeFastaFiles resources/genome_ncbi/GCF_000181335.3_Felis_catus_9.0_genomic.fna --sjdbGTFfile resources/genome_ncbi/GCF_000181335.3_Felis_catus_9.0_genomic.gtf --sjdbOverhang 149

## set up some relative paths
export DATA=data/trimmed
mkdir results/star_quant
export OUTDIR=results/star_quant

## Load genome index
STAR --genomeLoad LoadAndExit --genomeDir $GENOMEDIR

## Loop over all read files
cd $DATA

for i in $(ls *R1_001_trimmed.fq.gz | sed 's/\R1_001_trimmed.fq.gz//'); do echo STAR --runThreadN 22 --genomeDir $GENOMEDIR/GCF_000181335.3_Felis_catus_9.0_genomic --genomeLoad LoadAndKeep --readFilesIn ${i}R1_001_trimmed.fq.gz ${i}R2_001_trimmed.fq.gz --readFilesCommand zcat --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode GeneCounts --outFileNamePrefix $OUTDIR/$i;
done

## Remove genome index from memory
STAR --genomeLoad Remove --genomeDir $GENOMEDIR
