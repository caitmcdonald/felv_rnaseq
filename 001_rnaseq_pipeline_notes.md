RNAseq pipeline notes
================
Last Updated: 2021-08-11

  - [Overview](#overview)
      - [RNA-seq samples](#rna-seq-samples)
      - [The cat genome](#the-cat-genome)
  - [Pipeline](#pipeline)
      - [fastqc + multiqc](#fastqc-multiqc)
      - [Trimgalore](#trimgalore)
      - [STAR (`quantMode`)](#star-quantmode)
      - [OrthoFinder (RBH analysis)](#orthofinder-rbh-analysis)
      - [ID genes 1 Mb up and downstream of LTR integration
        sites](#id-genes-1-mb-up-and-downstream-of-ltr-integration-sites)
      - [edgeR](#edger)
      - [WGCNA](#wgcna)

## Overview

### RNA-seq samples

We have 38 RNA-seq libraries from the following samples:

| cat\_id  | cell\_type | population | status     | read | fastq                                        |
| :------- | :--------- | :--------- | :--------- | :--- | :------------------------------------------- |
| 4438     | PBMC       | SPF        | uninfected | R1   | 4438\_S1\_L002\_R1\_001.fastq.gz             |
| 4438     | PBMC       | SPF        | uninfected | R2   | 4438\_S1\_L002\_R2\_001.fastq.gz             |
| 4460     | PBMC       | SPF        | uninfected | R1   | 4460\_S2\_L002\_R1\_001.fastq.gz             |
| 4460     | PBMC       | SPF        | uninfected | R2   | 4460\_S2\_L002\_R2\_001.fastq.gz             |
| 4474     | PBMC       | SPF        | uninfected | R1   | 4474\_S3\_L002\_R1\_001.fastq.gz             |
| 4474     | PBMC       | SPF        | uninfected | R2   | 4474\_S3\_L002\_R2\_001.fastq.gz             |
| 4501     | PBMC       | SPF        | uninfected | R1   | 4501\_S4\_L002\_R1\_001.fastq.gz             |
| 4501     | PBMC       | SPF        | uninfected | R2   | 4501\_S4\_L002\_R2\_001.fastq.gz             |
| 4510     | PBMC       | SPF        | uninfected | R1   | 4510\_S5\_L002\_R1\_001.fastq.gz             |
| 4510     | PBMC       | SPF        | uninfected | R2   | 4510\_S5\_L002\_R2\_001.fastq.gz             |
| 4520     | PBMC       | SPF        | uninfected | R1   | 4520\_S6\_L002\_R1\_001.fastq.gz             |
| 4520     | PBMC       | SPF        | uninfected | R2   | 4520\_S6\_L002\_R2\_001.fastq.gz             |
| DC1      | fibroblast | outbred    | uninfected | R1   | DC1MINUS\_S7\_L002\_R1\_001.fastq.gz         |
| DC1      | fibroblast | outbred    | uninfected | R2   | DC1MINUS\_S7\_L002\_R2\_001.fastq.gz         |
| DC1      | fibroblast | outbred    | infected   | R1   | DC1PLUS\_S8\_L002\_R1\_001.fastq.gz          |
| DC1      | fibroblast | outbred    | infected   | R2   | DC1PLUS\_S8\_L002\_R2\_001.fastq.gz          |
| DC2      | fibroblast | outbred    | uninfected | R1   | DC2MINUES\_S9\_L002\_R1\_001.fastq.gz        |
| DC2      | fibroblast | outbred    | uninfected | R2   | DC2MINUES\_S9\_L002\_R2\_001.fastq.gz        |
| DC2      | fibroblast | outbred    | infected   | R1   | DC2PLUS\_S10\_L002\_R1\_001.fastq.gz         |
| DC2      | fibroblast | outbred    | infected   | R2   | DC2PLUS\_S10\_L002\_R2\_001.fastq.gz         |
| DC2      | fibroblast | outbred    | unknown    | I1   | DC2Pool\_NoIndex\_L001\_I1\_001.fastq.gz     |
| DC2      | fibroblast | outbred    | unknown    | R1   | DC2Pool\_NoIndex\_L001\_R1\_001.fastq.gz     |
| DC2      | fibroblast | outbred    | unknown    | R2   | DC2Pool\_NoIndex\_L001\_R2\_001.fastq.gz     |
| DC3      | fibroblast | outbred    | unknown    | I1   | DC3Pool\_NoIndex\_L001\_I1\_001.fastq.gz     |
| DC3      | fibroblast | outbred    | unknown    | R1   | DC3Pool\_NoIndex\_L001\_R1\_001.fastq.gz     |
| DC3      | fibroblast | outbred    | unknown    | R2   | DC3Pool\_NoIndex\_L001\_R2\_001.fastq.gz     |
| DC4      | fibroblast | outbred    | uninfected | R1   | DC4MINUS\_S11\_L002\_R1\_001.fastq.gz        |
| DC4      | fibroblast | outbred    | uninfected | R2   | DC4MINUS\_S11\_L002\_R2\_001.fastq.gz        |
| DC4      | fibroblast | outbred    | infected   | R1   | DC4PLUS\_S12\_L002\_R1\_001.fastq.gz         |
| DC4      | fibroblast | outbred    | infected   | R2   | DC4PLUS\_S12\_L002\_R2\_001.fastq.gz         |
| Mischief | fibroblast | puma       | uninfected | R1   | Mischief\_Minus\_S12\_L002\_R1\_001.fastq.gz |
| Mischief | fibroblast | puma       | uninfected | R2   | Mischief\_Minus\_S12\_L002\_R2\_001.fastq.gz |
| Mischief | fibroblast | puma       | infected   | R1   | Mischief\_Plus\_S11\_L002\_R1\_001.fastq.gz  |
| Mischief | fibroblast | puma       | infected   | R2   | Mischief\_Plus\_S11\_L002\_R2\_001.fastq.gz  |
| X2654    | fibroblast | outbred    | uninfected | R1   | X2654MINUS\_S13\_L002\_R1\_001.fastq.gz      |
| X2654    | fibroblast | outbred    | uninfected | R2   | X2654MINUS\_S13\_L002\_R2\_001.fastq.gz      |
| X2654    | fibroblast | outbred    | infected   | R1   | X2654PLUS\_S14\_L002\_R1\_001.fastq.gz       |
| X2654    | fibroblast | outbred    | infected   | R2   | X2654PLUS\_S14\_L002\_R2\_001.fastq.gz       |

Samples are named such as: `4438_S1_L002_R1_001.fastq.gz` in the
`data/raw` directory. L002 is common to all files and indicates
sequencing lane. 001 is also common to all files and is just the
standard Illumina append. S indicates sample. I could shorten all
filenames in `data` like so:

``` bash
for fname in *.fastq.gz ; do mv "$fname" "$(echo "$fname" | sed -r 's/L002_//')" ; done
for fname in *.fastq.gz ; do mv "$fname" "$(echo "$fname" | sed -r 's/_001//')" ; done
```

However, it’s best practice to **not** rename samples. For Snakemake, a
better option is to hardlink to the files, rename with tidy sample names
(e.g. S001), then use input functions and wildcards on these.

### The cat genome

The most recent cat genome assembly is:
[Felis\_catus\_9.0](https://www.ncbi.nlm.nih.gov/assembly/GCF_000181335.3/),
which has a total sequence length of 2.52 Gb and 19 chromosomes.

According to the [RefSeq Annotation
Report](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Felis_catus/104/),
it consists of:

  - Genes/psuedogenes: 35,588 (or else 31,417???)
      - 19,748 protein-coding
  - mRNAs: 54,713

<br>

## Pipeline

I ran the following pipeline with sbatch scripts or in R:

1.  run fastqc + multiqc
2.  trim
3.  re-run fastqc + multiqc on trimmed libraries
4.  run STAR in quantmode
5.  run OrthoFinder to restrict to cat/puma orthologs
6.  restrict to genes that fall within 1Mb of LTR integration sites
7.  edgeR
8.  WGCNA

### fastqc + multiqc

Because fastqc assigns 250MB per cpu and I want to try using 19 cores, I
requested 5GB (250 x 19 = 4750).

I submitted the
[fastqc\_multiqc\_slurm.sh](workflow/scripts/fastqc_multiqc_slurm.sh)
script via:

``` bash
sbatch workflow/scripts/fastqc_multiqc_slurm.sh
```

With the `--test-only` flag first to make sure I wasn’t requesting more
memory than existed. This took \~2.5 hours to run. This could definitely
be better optimized as a job array…

### Trimgalore

Trimgalore worked great as a loop with standard quality and read length
flags. See [trimgalore\_slurm.sh](scratch/trimgalore_slurm.sh).

### STAR (`quantMode`)

Ran STAR in `quantMode`, which generates gene-level quantification
files. [starquant\_slurm.sh](scratch/starquant_slurm.sh).

I then used bash one-liners to compile the **unstranded** output
(Elliott used an unstranded library prep) into an counts matrix
according to
[readcounts\_to\_matrix.sh](workflow/scripts/readcounts_to_matrix.sh)

### OrthoFinder (RBH analysis)

**Background:** If we want to compare gene expression between our cat
and panther samples, we need to restrict analysis to only orthologs.
There are a number of ways to do this. Classically, we do it with a
reciprocal best hit `blastp` run (blast all cat samples to puma genome,
blast all puma samples to cat genome, retain only genes with hits both
ways as presumptive orthologs). There are now many programs that will do
this for you that are much faster than a local blast alignment. One is
OrthoFinder, which uses `diamond`.

OrthoFinder requires protein fasta files for each proteome of interest,
then performs an all-vs-all blast with `diamond`. An added bonus is that
it also uses gene trees to infer orthology, which is more rigorous than
some other approaches.

I’m following this
[tutorial](https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html).

OrthoFinder creates a results dir called ‘OrthoFinder’ inside the input
proteome directory. For the RBH, we basically only care about
Orthologues output.

**Quality control:** Looking at the statistics files, we see that only
\~50% of cat genes are assigned to orthogroups whereas 90% of FP genes
are. Anything under 80% indicates poor species sampling. This
is…unexpected\! I would have assumed the opposite. Sue mentioned that
the cat genome was notoriously bad…(I asked Michelle and she said there
are lots of reasons for this: genotyping array is still only 100K when a
600K array has been promised for a decade, there’s not a lot of
money/industry money is in the hands of people who aren’t the best at
genetics, cats are all mutts so imputation is harder, etc.).

Other methods of QC are to look at the gene trees and the gene
duplication events to make sure they make sense…because we only have two
species, I don’t think the trees are really going to show anything
useful.

**Ortholog file:** Since most of our samples are cats and we only have
one panther, it makes more sense to use the cat ortholog file than the
panther ortholog file. I can filter the respective read counts matrices
to only orthologs, and combine for DGE analysis.

The ortholog file lists protein IDs, and I need to convert to gene IDs:

The easiest way to do this is to read the gtf file into R and use
`rtracklayer`:

``` r
library(rtracklayer)
gtf <- rtracklayer::import("results/orthofinder/GCF_000181335.3_Felis_catus_9.0_genomic.gtf")
gene_tran_pro <- as.data.frame(mcols(gtf)[,c("gene_id","transcript_id","protein_id")])
```

Then, I can get a list of cat orthologs (protein\_id) from Orthofinder
output (.csv file) and join to locate gene\_id:

``` r
orthologs <- read_delim("results/orthofinder/GCF_000181335.3_Felis_catus_9.0_protein__v__GCF_003327715.1_PumCon1.0_protein.csv", delim="\t", col_names = c("Orthogroup","Cat_protein_id","puma_protein_id"), skip = 1) %>% 
  separate_rows(Cat_protein_id, convert=TRUE, sep = ", ")

ortho_and_gene <- left_join(orthologs, gene_tran_pro, by=c("Cat_protein_id"="protein_id")) %>%
  select(gene_id) %>% 
  distinct()
```

There are multiple proteins derived from each gene, so we need to drop
repeats

``` r
dat <- read_delim("readcountsmatrix.txt", delim = "\t") %>% select(-DC2Pool, -DC3Pool)
counts <- dat[-c(1:4, 31498),]

counts_w_protein <- left_join(ortho_and_gene, counts, by=c("gene_id"="genes"))

counts_ortho_only <- counts_w_protein %>% 
  #select(-Cat_protein_id, -transcript_id) %>% 
  distinct() %>% 
  drop_na()
# write_delim(counts_ortho_only, "results/orthofinder/genecounts_orthologs.txt", delim="\t")
```

### ID genes 1 Mb up and downstream of LTR integration sites

After restricting to orthologs, we can further restrict our gene set to
only those genes that are within 1 Mb of the \~700 LTR integration sites
that Elliott identified.

Using [Supplemental Table
1](data/elliot_ltr_manuscript/Supplemental%20Table%201%20Run%202%20Curated%20integration%20sites.xlsx)
from Elliott’s integration site manuscript, I extracted LTR start sites,
then [generated a list of chromosomal regions +/- 1
Mb](data/LTR_allsites_upanddown.txt).

I used this list of defined regions (which is in bed format) to query
the [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables) to
extract all genes within these windows in bed format.

``` r
LTR_genes <- read_delim("data/felCat9_1Mb_LTRsites_all.txt", delim="\t", skip = 1, col_names = F) %>% 
  select(X4) %>% 
  distinct()

LTR_genes %>% summarise(n())
```

Based on the genes extracted from UCSC, it appears there are 31,076
genes within 1 Mb of an LTR integration site. This is essentially the
entire genome\! Thus, it’s not going to give us any more specific
results than we already have. Options for parsing this out better:

**1. Restrict to \<1 Mb up and downstream.**

    - The 1 Mb decision was because this is the maximum distance for LTR enhancer function. We could instead look just at promoter function? What would this distance be?

**2. Pick out particularly interesting integration sites and look only
at those.**

    - For example, only the integration sites found in >10 cats
    - __I'll use this second approach__

#### Interesting integration sites

**1. Fibroblasts: infected vs. uninfected** Restrict to LTR sites
present in ≥3 (out of 4) fibroblast samples = [78 LTR
sites](data/LTR_fibro_upanddown.txt)

``` r
gtf <- rtracklayer::import("results/orthofinder/GCF_000181335.3_Felis_catus_9.0_genomic.gtf")
gene_tran_pro <- as.data.frame(mcols(gtf)[,c("gene_id","transcript_id","protein_id")])

LTR_fibro <- read_delim("data/felCat9_1Mb_LTRsites_fibro.txt", delim="\t", skip = 1, col_names = F) %>% 
  select(X4) %>% 
  distinct()
LTR_fibro_genid <- left_join(LTR_fibro, gene_tran_pro, by=c("X4"="transcript_id")) %>% 
  select(-protein_id, -X4) %>% 
  distinct()

# write_tsv(LTR_fibro_genid, "data/LTR_fibro_genids.txt")
```

This doesn’t yield any interesting genes…

**2. PMBCs: present vs absent** Restrict to LTR sites present in 3 (out
of 6) PMBC samples = [42 sites](data/LTR_pmbc_upanddown.txt). These need
to be partitioned into three separate groups for comparison in edgeR.

When we run a DGE analysis on these partitioned groups, we get:

    __Set1:__ 15 sig genes; no sig LTR genes
    __Set2:__ 51 sig genes; not sig LTR genes
    __Set3:__ 1 significant gene; no significant LTR genes

**3. Present in ≥10 cats** Restrict to sites only found in ≥10 cats in
Elliott’s larger study = [87 sites](data/LTR_over10cats_upanddown.txt)

But then for this, what groups would I compare? There are no clear
groupings where this would work, except for PMBCs vs. fibroblasts, which
is super confounded…

**3a. Compare PMBCs vs fibroblasts:** When restricting to genes 1Mb from
LTR insertion sites present in \>10 cats, there are 687 significant
genes. But, this isn’t meaningful at all\! The sites that are present in
\>10 cats may or may not be present in the PMBC and fibroblast
samples.If they aren’t present in one or the other, then the expression
could be due to cell type or LTR presence/absence and there’s no way to
tell. I could try to: a) restrict to only LTR sites present in BOTH, b)
bin PMBC+fibro present vs PMBC+fibro absent and re-run. However, there
are not enough samples to yield at least 3 bio reps per treatment in
this design, so it’s also not an option.

### edgeR

I used the [readcountsmatrix.txt](readcountsmatrix.txt) and the sample
[metadata](data/felv_metadata.tsv) to begin edgeR analyses.

**Note:** when restricting to smaller gene sets (e.g. 1-3 above), all
normalization and model fitting must be done on the full gene set. Then
at the T-test to ID significant genes, I can restrict to genes of
interest to recalculate the P-values. So, because there are no
significantly differentially expressed genes between infected and
uninfected fibroblasts, there will be no significant genes from the
restricted set \#1…

#### Data exploration

Code: [DGE\_dataexplore.R](workflow/scripts/DGE_dataexplore.R)

Looking at all samples, we can see that there are 31,498 genes in the
full dataset. After modest filtering for sequence errors, that goes down
to \~16,000 genes. Library size ranged from 2,692,910 (MischiefPlus) to
33,912,615 (DC4PLUS). That’s way too few reads for the Mischief (puma)
samples. When we look at sample correlation and a PCA, we can also see
that they’re way different from the other fibroblast samples, which is
probably a function of poor sequencing, and is not biologically
relevant. This relationship remains, even when we filter to only
cat-puma orthologs. Consider removing the Mischief samples from further
analyses.

#### Differential expression

Code: [DGE\_edgeR.R](workflow/scripts/DGE_edgeR.R)

So there are a couple comparisons we can make:

  - Cell type differential expression (PMBCs vs fibroblasts)
  - uninfected vs. FeLV-infected differential expression (fibroblasts
    only)
  - LTR-specific:
      - LTR+ vs LTR- PMBCs
      - LTR+ vs LTR- fibroblasts, inf vs uninf

Looking at cell type, 4510 is a big outlier, again I suspect due to poor
sequencing, although Elliott says this animal has far more LTR sites
than any of the other samples. We can’t really disentangle this because
it did have far fewer reads sequenced than the other individuals, so it
could either be sequencing or LTR number.

Expression is clearly cell-specific, which we would expect. There are
\~6,000 genes significantly differentially expressed between fibroblasts
and PMBCs.

Looking at FeLV infection status, we see that infection status plays a
pretty negligible role. Samples cluster according to biological
replicate and not according to infection status. Thus, we end up with
only 3 significantly DE genes, all of which are uncharacterized
proteins.

### WGCNA

We could try running WGCNA to look at how LTR copy number affects
expression because WGCNA is good for correlating expression with
continuous variables. However, we can’t use WGCNA for sample sizes \<15…
