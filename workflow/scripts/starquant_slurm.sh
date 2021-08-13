#!/bin/bash

#SBATCH --mem 110GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=coby.mcdonald@colostate.edu
#SBATCH --output=logs/star/output-%j
#SBATCH --error=logs/star/error-%j

export WORK=/scratch/summit/camcd@colostate.edu/felv_rnaseq

## Get genome and gtf
# export NCBI=$WORK/resources/genome_ncbi
# cd $NCBI
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_000181335.3_Felis_catus_9.0/GCF_000181335.3_Felis_catus_9.0_genomic.fna.gz
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Felis_catus/latest_assembly_versions/GCF_000181335.3_Felis_catus_9.0/GCF_000181335.3_Felis_catus_9.0_genomic.gtf.gz
# gunzip *.gz


## Set relative paths
# mkdir $WORK/resources/genome_star
export STARDIR=$WORK/resources/genome_star
export FASTA=$WORK/resources/genome_ncbi/GCF_000181335.3_Felis_catus_9.0_genomic.fna
export GTF=$WORK/resources/genome_ncbi/GCF_000181335.3_Felis_catus_9.0_genomic.gtf
export SAMPLES=$WORK/data/trimmed
export OUTDIR=$WORK/results/star_quant_bam


## Run STAR with quant mode
cd $WORK
source ~/.bashrc
conda activate star_quant

# Build genome index
STAR --runThreadN 22 --runMode genomeGenerate --genomeDir $STARDIR --genomeFastaFiles $FASTA --sjdbGTFfile $GTF --sjdbOverhang 149

# Load genome index
STAR --genomeLoad LoadAndExit --genomeDir $STARDIR

## Loop over all read files
cd $SAMPLES

for i in $(ls *R1_001_val_1.fq.gz | sed 's/\R1_001_val_1.fq.gz//'); do STAR --runThreadN 22 --genomeDir $STARDIR --readFilesIn ${i}R1_001_val_1.fq.gz ${i}R2_001_val_2.fq.gz --readFilesCommand zcat --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix $OUTDIR/$i;
done

#--quantMode TranscriptomeSAM GeneCounts will return aligned bams and gene counts

## Remove genome index from memory
STAR --genomeLoad Remove --genomeDir $STARDIR
