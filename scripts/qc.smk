import glob, sys #following Eric A--omg PYTHON exciting!!!

fullnames= glob.glob('/data/*.fastq.gz')

# get sample name without .fastq.gz
sample= []
for name in fullnames:
    base= name.split('.')[0]
    sample.append(base)

# get base name without R1 or R2
basename= []
for read in sample:
    noR= read.split('_')[0] + '_' + read.split('_')[1]
    basename.append(noR)

# don't totally understand this rule yet...
#rule all:
#    input:
        #the last output

rule fastqc_raw:
    input:
        "/data/{sample}.fastq.gz"
    output:
        "/results/fastqc/raw/{sample}.html",
        "/results/fastqc/raw/{sample}.zip"
    conda:
        "felv_rna.yml"
    log:
        "logs/fastqc/raw/{sample}.log"
    shell:
        "fastqc {input}"

rule multiqc_raw:
    input:
        "/results/fastqc/raw/"
    output:
        "/results/multiqc/raw/multiqc_report.html",
        #directory("/results/multiqc/raw/multiqc_data")
    shell:
        "multiqc {input}"

rule trimgalore:
    input:
        "/data/{sample}.fastq.gz"
    output:
        "/results/trimming/{basename}_R1.fastq.gz_trimming_report.txt",
        "/results/trimming/{basename}_val_1.fastq.gz",
        "/results/trimming/{basename}_R2.fastq.gz_trimming_report.txt",
        "/results/trimming/{basename}_val_2.fastq.gz"
    shell:
        "trim_galore -q 5 --paired --cores 2 {input}"

rule fastqc_trim:
    input:
        "/results/trimming/{sample}.trimmed.fastq.gz"
    output:
        "results/fastqc/trimmed/{sample}.trimmed_fastqc.zip",
        "results/fastqc/trimmed/{sample}.trimmed_fastqc.html"

rule multiqc_trim:
    input:
        "/results/fastqc/trimmed/"
    output:
        "/results/multiqc/trimmed/multiqc_report.html",
        #directory("/results/multiqc/trimmed/multiqc_data")
