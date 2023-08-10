

wd="/home/amovas/data/genome-evo-proj"
EXP,LINES,SAMPLES = glob_wildcards("/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz")



# If you want to run the pipeline for a subset of samples use these assignments:

#EXP = ["iv","iv", "iii", "iii"]
#LINES = ["20","20", "15", "15"]
#SAMPLES = ["20MT2EXPIVVP180seq09052019_S24_L001","20MT2EXPIVVP410seq08122022_S33_L001", "15MT4EXPIIIVP440seq20052022_S16_L001", "15MT4EXPIIIVP510seq09092022_S9_L001"]





rule all:
    input:
       "/home/amovas/data/genome-evo-proj/results/tables/3-p/all_variants.pkl",
       expand("/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}_aligned.bam.bai", zip,experiment=EXP, line=LINES, sample=SAMPLES)



rule bwa_map:
    input:
        "/home/amovas/data/genome-evo-proj/data/reference/plasmid/hiv_plasmid_ref_genome.fasta", #ancestor/ancestor_consensus.fasta
        "/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz",
        "/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R2_001.fastq.gz"
    output:
        temp("/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}.bam.temp")
    log:
        "/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/logs/{experiment}_{line}_{sample}.log"
    threads: 1
    shell:
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"



rule samtools_sort:
    input:
        "/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}.bam.temp"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}_aligned.bam"
    shell:
        "samtools sort {input} > {output}"



rule samtools_index:
    input:
        "/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}_aligned.bam"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}_aligned.bam.bai"
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
        bam="/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}_aligned.bam",
        bai="/home/amovas/data/genome-evo-proj/data/processed-data/mappings/3-p/{experiment}/{line}/{sample}_aligned.bam.bai"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/mutations/3-p/{experiment}/{line}/{sample}.csv"
    shell:
        "python /home/amovas/data/genome-evo-proj/src/snakemake-run/my_call_mutations.py {input.bam} /home/amovas/data/genome-evo-proj/data/reference/plasmid/hiv_plasmid_ref_genome.fasta AF324493.2 > {output}"

# AF324493.2
# NL43seq210314
# K03455.1
# NL43_ann_wk0virusPassRef_SEQregion



rule collect_mutations:
    input:
        expand("/home/amovas/data/genome-evo-proj/data/processed-data/mutations/3-p/{experiment}/{line}/{sample}.csv",
                zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        #"mutations/all_ins.pkl",
        #"mutations/all_dels.pkl",
        "/home/amovas/data/genome-evo-proj/results/tables/3-p/all_variants.pkl"
    shell:
        "python /home/amovas/data/genome-evo-proj/src/snakemake-run/my_collect_mutations.py {input}"

