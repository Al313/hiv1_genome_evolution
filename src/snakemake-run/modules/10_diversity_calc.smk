

rule calculateDiversity:
    input:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_annotated_variants.tsv.gz"
    output:
        f"{wd}/results/tables/misc/diversity/{{line}}_diversity_measures.tsv
    shell:
        "Rscript {wd}/src/analysis/diversity-calculations/get_diversity.R {input}"
