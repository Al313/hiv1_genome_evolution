rule prepareFeatureList:
    input:
        f"{wd}/data/reference/annotations/features/sequence-features.tsv"
    output:
        f"{wd}/results/tables/pipeline-outputs/cds_feature_list.tsv"
    script:
        f"{wd}/src/snakemake-run/R-scripts/prep_feature_list.R"


rule extractAnnotation:
    input:
        f"{wd}/results/tables/pipeline-outputs/cds_feature_list.tsv",
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.ann.vcf.gz"
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_annotated_variants.tsv.gz"
    resources:
        mem_mb=50000
    params:
        wd = wd
    shell:
        "Rscript {wd}/src/snakemake-run/R-scripts/aa_change.R {wildcards.experiment} && Rscript {wd}/src/snakemake-run/R-scripts/add_mut_context.R {wildcards.experiment}"

