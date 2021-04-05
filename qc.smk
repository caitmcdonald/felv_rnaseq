# import glob, sys #following Eric A--omg PYTHON exciting!!!
#
# fullnames= glob.glob('/data/*.fastq.gz')
#
# # get sample name without .fastq.gz
# sample= []
# for name in fullnames:
#     base= name.split('.')[0]
#     sample.append(base)
#
# # get base name without R1 or R2
# basename= []
# for read in sample:
#     noR= read.split('_')[0] + '_' + read.split('_')[1]
#     basename.append(noR)

# Import metadata
# import pandas as pd
#
# metadata = pd.read_table("felv_metadata_long.txt").set_index("sample_id", drop=False)

configfile: "config.yaml"

# don't totally understand this rule yet...
rule all:
    input:
        #the last output
        "/results/multiqc/trimmed/multiqc_report.html"

rule fastqc_raw:
    input:
        # "/data/{sample}.fastq.gz"
        # get_fastqc_input_fastqs
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "/results/fastqc/raw/{sample}.html",
        "/results/fastqc/raw/{sample}.zip"
    conda:
        "qctrim.yml"
    log:
        "logs/fastqc/raw/{sample}.log"
    threads: 1
    shell:
        "(fastqc {input}) 2> {log}"

# rule multiqc_raw:
#     input:
#         "/results/fastqc/raw/"
#     output:
#         "/results/multiqc/raw/multiqc_report.html",
#         #directory("/results/multiqc/raw/multiqc_data")
#     shell:
#         "multiqc {input}"

# rule cutadapt:
#     input:
#         "/data/{sample}.fastq.gz"
#     output:
#     shell:
#
# rule fastqc_trim:
#     input:
#         "/results/trimming/{basename}_val_1.fastq.gz",
#         "/results/trimming/{basename}_val_2.fastq.gz"
#     output:
#         "/results/fastqc/trimmed/{sample}.trimmed_fastqc.zip",
#         "/results/fastqc/trimmed/{sample}.trimmed_fastqc.html"
#     shell:
#         "fastqc {input}"
#
# rule multiqc_trim:
#     input:
#         "/results/fastqc/trimmed/"
#     output:
#         "/results/multiqc/trimmed/multiqc_report.html",
#         #directory("/results/multiqc/trimmed/multiqc_data")
#     shell:
#         "multiqc {input}"
