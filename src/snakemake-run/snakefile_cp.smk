


wd="/Users/alimos313/Documents/studies/phd/research/genome-evo-proj"
EXP,LINES,SAMPLES = glob_wildcards("/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz")

# run the pipeline for a subset of samples:
subset_run = True
if subset_run:
    EXP=EXP[::52]
    LINES=LINES[::52]
    SAMPLES=SAMPLES[::52]
    final_out = "all_muts_subset.pkl"
else:
    final_out = "all_muts.pkl"


rule all:
    input:
       f"/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mutations/{final_out}",
#       expand("/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}.bam.temp", zip,experiment=EXP, line=LINES, sample=SAMPLES)
#       expand("/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}_aligned.bam.bai", zip,experiment=EXP, line=LINES, sample=SAMPLES)



rule bwa_map:
    input:
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/reference2/HIV-ref-genome.fasta", #ancestor/ancestor_consensus.fasta
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz",
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R2_001.fastq.gz"
    output:
        temp("/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}.bam.temp")
    log:
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/logs/{experiment}_{line}_{sample}.log"
    threads: 1
    shell:
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"



rule samtools_sort:
    input:
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}.bam.temp"
    output:
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}_aligned.bam"
    shell:
        "samtools sort {input} > {output}"



rule samtools_index:
    input:
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}_aligned.bam"
    output:
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}_aligned.bam.bai"
    shell:
        "samtools index {input}"



"""
rule quality_assessment:
    input:
        bam="HIV_LTE_mappings/{experiment}/rep/{sample}_aligned.bam",
        bai="HIV_LTE_mappings/{experiment}/rep/{sample}_aligned.bam.bai"
    output:
        "quality_control/coverage_plots/{experiment}-{line}-{sample}.png",
        "quality_control/{experiment}/rep/{sample}.quals"
    script:
        "snakemake_scripts/quality_control_mapping.py"

rule collect_quality_assessment:
    input:
        #expand("quality_control/{exp}/{line}/{sample}.quals",
        #       zip, exp=directories,line=lines, sample=IDS),
        expand("quality_control/{experiment}/rep/{sample}.quals",
                zip, line=lines2, sample=IDS2)
    output:
        "quality_control/all.quals"
    shell:
        "cat {input} >> {output}"

"""



rule bam_to_mutations:
    input:
        bam="/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}_aligned.bam",
        bai="/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mapped/{experiment}/{line}/{sample}_aligned.bam.bai"
    output:
        "/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mutations/{experiment}/{line}/{sample}.csv"
    shell:
        "python /Users/alimos313/Documents/studies/phd/research/genome-evo-proj/code/snakemake-scripts/call_mutations.py {input.bam} /Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/reference2/HIV-ref-genome.fasta AF324493.2 > {output}"

# AF324493.2
# NL43seq210314
# K03455.1
# NL43_ann_wk0virusPassRef_SEQregion



rule collect_mutations:
    input:
        expand("/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mutations/{experiment}/{line}/{sample}.csv",
                zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        #"mutations/all_ins.pkl",
        #"mutations/all_dels.pkl",
        f"/Users/alimos313/Documents/studies/phd/research/genome-evo-proj/data/processed-data/mutations/{final_out}"
    shell:
        "echo {input} & python /Users/alimos313/Documents/studies/phd/research/genome-evo-proj/code/snakemake-scripts/collect_mutations.py {input}"

