
rule getMappingQual:
    input:
        bam=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    output:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted_stats/qualimapReport.html"
    shell:
        "qualimap bamqc -bam {input.bam}"


rule getCoverage:
    input:
        bam=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam",
        bai=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam.bai"
    output:
        f"{wd}/data/processed-data/coverage-assessment/pipeline-outputs/coverage-plots/{{experiment}}-{{line}}-{{sample}}.png",
        f"{wd}/data/processed-data/coverage-assessment/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.quals"
    script:
        "{wd}/src/snakemake-run/python-scripts/quality_control_mapping.py"

rule collectCoverage:
    input:
        #expand("quality_control/{exp}/{{line}}/{{sample}}.quals",
        #       zip, exp=directories,line=lines, sample=IDS),
        expand(f"{wd}/data/processed-data/coverage-assessment/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.quals",
                zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        f"{wd}/results/tables/pipeline-outputs/qc/all_quals.tsv.gz"
    shell:
        "cat {input} | gzip >> {output}"

