rule generateVCF:
    input:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.csv.gz"
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.vcf.gz"
    shell:
        "Rscript {wd}/src/pipeline/R-scripts/variant_to_vcf.R {wildcards.experiment} && gzip {wd}/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.vcf"



rule annotateVCF:
    input:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.vcf.gz"
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.ann.vcf.gz"
    shell:
        "cd {wd}/data/reference/annotations && snpEff -v hiv_plasmid_ref_genome {wd}/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.vcf.gz > {wd}/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.ann.vcf \
	&& gzip {wd}/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.ann.vcf"


