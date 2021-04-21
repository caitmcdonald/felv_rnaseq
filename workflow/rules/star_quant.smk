## Run STAR with quant mode

rule star_index:
    input:
    output:
    conda:
        "envs/star_quant.yaml"
    log:
    shell:

"echo STAR --runThreadN 22 --runMode genome Generate --genomeDir {input.index} genomeFastaFiles {input.fasta} --sjdbGTFfile {input.gtf} --sjdbOverhang 149"

rule star_quant:
    input:
    output:
    conda:
        "envs/star_quant.yaml"
    log:
    shell:

"echo STAR --runThreadN 22 --genomeDir --readFilesIn {input.fq1} {input.fq2}  --readFilesCommand zcat --outFilterMismatchNmax 3 --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterIntronMotifs RemoveNoncanonicalUnannotated --quantMode GeneCounts"
