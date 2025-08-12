rule callMutation:
    input:
        bam=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam",
        bai=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam.bai"
    output:
        f"{wd}/data/processed-data/mutations/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.csv"
    shell:
        "python {wd}/src/snakemake-run/python-scripts/call_mutations.py {input.bam} {wd}/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta NL43_ann_wk0virusPassRef_plasmid > {output}"


# Helper function to get samples for a specific experiment
def get_samples_for_experiment(wildcards):
    exp_samples = []
    for i, exp in enumerate(EXP):
        if exp == wildcards.experiment:
            exp_samples.append(f"{wd}/data/processed-data/mutations/pipeline-outputs/{EXP[i]}/{LINES[i]}/{SAMPLES[i]}.csv")
    return exp_samples


rule collectMutation:
    input:
        get_samples_for_experiment
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.csv.gz"
    shell:
        "python {wd}/src/snakemake-run/python-scripts/collect_mutations.py {input}"


