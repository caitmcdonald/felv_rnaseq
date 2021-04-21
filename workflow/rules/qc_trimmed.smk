## Run fastqc on trimmed files, then run multiqc on fastqc output
import glob, sys
trimmed = glob.glob('data/trimmed/*fastq.gz')

filename = []
for f in trimmed:
  no_ext = f.split('.')[0]
  filename.append(no_ext)

rule fastqc_trimmed:
    input:
        trimmed
    output:
        "results/fastqc/trimmed/{filename}_fastqc.html",
        "results/fastqc/trimmed/{filename}_fastqc.zip"
    conda:
        "envs/fastqc.yaml" #note: don't list workflow parent dir as snakemake recognizes it automatically
    log:
        "logs/fastqc/trimmed/{filename}.log"
    threads: 19
    shell:
        "echo (fastqc -t {threads} {input}) 2> {log}"

rule multiqc_trimmed:
    input:
        "results/fastqc/trimmed"
    output:
        "results/multiqc/trimmed/multiqc_report.html",
        direct=directory("results/multiqc/trimmed")
    conda:
        "envs/multiqc.yaml"
    shell:
        "echo multiqc {input} -o {output.direct}"
