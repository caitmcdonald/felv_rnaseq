## Run fastqc on raw files, then run multiqc on fastqc output
import glob, sys
raw = glob.glob('data/raw/*fastq.gz')

filename = []
for f in raw:
    no_ext = f.split('.')[0]
    filename.append(no_ext)

fastqc_raw_out = []
for raw in filename:
    new_filename = 'results/fastqc/raw/' + raw.split('/')[2] + '_fastqc.html'
    fastqc_raw_out.append(new_filename)

rule fastqc_raw:
    input:
        raw
    output:
        "results/fastqc/raw/{filename}_fastqc.html",
        "results/fastqc/raw/{filename}_fastqc.zip"
    conda:
        "../envs/fastqc.yaml" #note: don't list workflow parent dir as snakemake recognizes it automatically
    log:
        "logs/fastqc/raw/{filename}.log"
    threads: 20
    shell:
        "(fastqc -t {threads} {input}) 2> {log}"

rule multiqc_raw:
    input:
        # "results/fastqc/raw" #for some reason can't use dir as input???
        fastqc_raw_out
    output:
        "results/multiqc/raw/multiqc_report.html",
        direct=directory("results/multiqc/raw")
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc {input} -o {output.direct}"
