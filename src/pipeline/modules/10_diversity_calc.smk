

rule calculateDiversity:
    input:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_annotated_variants.tsv.gz"
    output:
        f"{wd}/results/tables/diversity/{{experiment}}/{{line}}_diversity_measures.tsv"
    shell:
        "Rscript /home/amovas/data/genome-evo-proj/src/pipeline/R-scripts/get_diversity.R {wildcards.line}"
