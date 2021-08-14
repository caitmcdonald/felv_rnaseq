enFeLV-exFeLV RNAseq: data and code repo
================
Coby McDonald

**last updated:** 2021-08-13

# Overview

## Project goal

Following work by [Chiu et
al. 2020](https://journals.asm.org/doi/full/10.1128/JVI.01274-20),
endogenous FeLV can potentially decrease the effects of exogenous FeLV
infection in domestic cats. The goal of this project is to identify
whether enFeLV LTR integration site significantly alters transcription
of host genes.

## Contents

This repo will include code, data, and results related to RNA-sequencing
of fibroblasts and PMBCs from domestic cats and one puma. PMBCs were
never infected with FeLV, so are included for baseline data. Fibroblast
samples are matched and include infected and uninfected samples. Running
analysis notes can be found in [RNAseq pipeline
notes](001_rnaseq_pipeline_notes.md).

## Pipeline

1.  Quality control with [fastQC +
    multiQC](workflow/scripts/fastqc_multiqc_raw_slurm.sh)
2.  Trimming with [trimgalore](workflow/scripts/trimgalore_slurm.sh)
3.  Repeat QC
4.  Quasi-alignment and quantification via
    [STAR](workflow/scripts/starquant_slurm.sh)
5.  RBH analysis with
    [OrthoFinder](workflow/scripts/orthofinder_slurm.sh)
6.  Differential gene expression
    [exploration](workflow/scripts/DGE_dataexplore.R) and
    [analysis](workflow/scripts/DGE_edgeR.R)
7.  Testing [WGCNA](workflow/scripts/wgcna.R) (but aborted this).
