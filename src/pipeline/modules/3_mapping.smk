rule bwaMapping:
    input:
        f"{wd}/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta", #ancestor/ancestor_consensus.fasta
        f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P.fastq.gz",
        f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P.fastq.gz"
    output:
        temp(f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.bam.temp")
    log:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/logs/{{experiment}}_{{line}}_{{sample}}.log"
    threads: 1
    shell:
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"



rule samtoolsSorting:
    input:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.bam.temp"
    output:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    shell:
        "samtools sort {input} > {output}"




rule samtoolsIndexing:
    input:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    output:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam.bai"
    shell:
        "samtools index {input}"

