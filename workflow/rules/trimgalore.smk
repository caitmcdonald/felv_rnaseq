## Trim with trimgalore for PE mode
rule trimgalore_pe:
    input:
        fq1=get_raw_fq1,
        fq2=get_raw_fq2
    output:
        # direct=directory("trimmed"),
        "data/trimmed/{sample_id}_R1_001_val_1.fastq.gz",
        "data/trimmed/{sample_id}_R2_001_val_2.fastq.gz",
        "data/trimmed/{sample_id}_R1_001.fastq.gz_trimming_report.txt",
        "data/trimmed/{sample_id}_R2_001.fastq.gz_trimming_report.txt",
        "data/trimmed/{sample_id}_R1_001_unpaired_1.fastq.gz",
        "data/trimmed/{sample_id}_R2_001_unpaired_2.fastq.gz"
    conda:
        "../envs/trim_galore.yaml"
    log:
        "logs/trimgalore/{sample_id}.log"
    threads: 20
    shell:
        "echo trim_galore --cores 20 --paired --retain_unpaired --phred33 --length 36 -q 5 --stringency 1 -e 0.1 -o /data/trimmed/ {input.fq1} {input.fq2}"
