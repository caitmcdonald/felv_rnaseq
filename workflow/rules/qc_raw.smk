## Run fastqc on raw files, then run multiqc on fastqc output
import glob, sys
raw = glob.glob('data/raw/*fastq.gz')

filename = []
for f in raw:
  no_ext = f.split('.')[0]
  filename.append(no_ext)

rule fastqc_raw:
    input:
        get_raw_fastqs
    output:
        "results/fastqc/raw/{filename}_fastqc.html",
        "results/fastqc/raw/{filename}_fastqc.zip"
    conda:
        "envs/fastqc.yaml" #note: don't list workflow parent dir as snakemake recognizes it automatically
    log:
        "logs/fastqc/raw/{filename}.log"
    threads: 19
    shell:
        "echo (fastqc -t {threads} {input}) 2> {log}"

rule multiqc_raw:
    input:
        "results/fastqc/raw"
    output:
        "results/multiqc/raw/multiqc_report.html",
        direct=directory("results/multiqc/raw")
    conda:
        "envs/multiqc.yaml"
    shell:
        "echo multiqc {input} -o {output.direct}"