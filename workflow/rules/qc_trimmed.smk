## Run fastqc on trimmed files, then run multiqc on fastqc output
import glob, sys
trimmed = glob.glob('trimmed/*val*')

filename = []
for file in trimmed:
  no_ext = file.split('.')[0]
  filename.append(no_ext)

fastqc_trim_out = []
for trim in filename:
  new_filename = 'results/fastqc/trimmed/' + trim.split('/')[2] + '_fastqc.html'
  fastqc_trim_out.append(new_filename)

rule fastqc_trimmed:
    input:
        trimmed
    output:
        "results/fastqc/trimmed/{filename}_fastqc.html",
        "results/fastqc/trimmed/{filename}_fastqc.zip"
    conda:
        "../envs/fastqc.yaml" #note: don't list workflow parent dir as snakemake recognizes it automatically
    log:
        "logs/fastqc/trimmed/{filename}.log"
    threads: 20
    shell:
        "(fastqc -t {threads} {input}) 2> {log}"

rule multiqc_trimmed:
    input:
        # "results/fastqc/trimmed"
        fastqc_trim_out
    output:
        "results/multiqc/trimmed/multiqc_report.html",
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc {input} -o results/multiqc/trimmed"
