rule getSeqQual:
    input:
        read1=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R1_001.fastq.gz",
        read2=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R2_001.fastq.gz"
    output:
        f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_R1_001_fastqc.html",
        f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_R2_001_fastqc.html"
    shell:
        "fastqc {input.read1} -o {wd}/data/processed-data/fastqc-reports/pipeline-outputs/{wildcards.experiment}/{wildcards.line} && fastqc {input.read2} -o {wd}/data/processed-data/fastqc-reports/pipeline-outputs/{wildcards.experiment}/{wildcards.line}"

# Helper function to get all fastqc outputs for a specific experiment and line
def get_fastqc_outputs_for_experiment_line(wildcards):
    fastqc_outputs = []
    for i, (exp, line) in enumerate(zip(EXP, LINES)):
        if exp == wildcards.experiment and line == wildcards.line:
            sample = SAMPLES[i]
            fastqc_outputs.extend([
                f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{exp}/{line}/{sample}_R1_001_fastqc.html",
                f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{exp}/{line}/{sample}_R2_001_fastqc.html"
            ])  
    return fastqc_outputs

rule collectSeqQual:
    input:
        get_fastqc_outputs_for_experiment_line
    output:
        f"{wd}/results/tables/pipeline-outputs/qc/{{experiment}}_{{line}}_multiqc_report.html"
    params:
        input_dir=f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{{experiment}}/{{line}}"
    shell:
        "cd {params.input_dir} && multiqc . -o multiqc && cp multiqc/multiqc_report.html {output} && rm -rf multiqc"


