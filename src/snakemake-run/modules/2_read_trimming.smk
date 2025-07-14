rule trimRead:
    input:
        one=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R1_001.fastq.gz",
        two=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R2_001.fastq.gz"
    output:
        one=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P",
        two=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P"
    shell:
        "x={output.one} && y=${{x%%_1P*}} && trimmomatic PE -phred33 {input.one} {input.two} -baseout ${{y}} LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36"

rule trimReadPostProcessing:
    input:
        one=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P",
        two=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P"
    output:
        one=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P.fastq.gz",
        two=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P.fastq.gz"
    shell:
        "gzip -S '.fastq.gz' {input.one} && gzip -S '.fastq.gz' {input.two}"

