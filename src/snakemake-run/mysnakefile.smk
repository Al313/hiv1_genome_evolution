

wd="/home/amovas/data/genome-evo-proj"
EXP,LINES,SAMPLES = glob_wildcards("/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz")



# If you want to run the pipeline for a subset of samples use these assignments:

#EXP = ["iv","iv", "iii", "iii"]
#LINES = ["20","20", "15", "15"]
#SAMPLES = ["20MT2EXPIVVP180seq09052019_S24_L001","20MT2EXPIVVP410combinedseq08122022_S33_L001", "15MT4EXPIIIVP440seq20052022_S16_L001", "15MT4EXPIIIVP510seq09092022_S9_L001"]





rule all:
    input:
       "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_annotated_variants.tsv.gz",
       "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/qc/all_quals.tsv.gz",
       expand("/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/qc/{experiment}_{line}_multiqc_report.html", zip, experiment=EXP, line=LINES),
       expand("/home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{experiment}/{line}/{sample}_R1_001_fastqc.html", zip,experiment=EXP, line=LINES, sample=SAMPLES)


rule sequencing_quality:
    input:
        read1="/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz",
        read2="/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R2_001.fastq.gz"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{experiment}/{line}/{sample}_R1_001_fastqc.html",
        "/home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{experiment}/{line}/{sample}_R2_001_fastqc.html"
    shell:
        "fastqc {input.read1} -o /home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{wildcards.experiment}/{wildcards.line} && fastqc {input.read2} -o /home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{wildcards.experiment}/{wildcards.line}"

rule collect_sequencing_quality:
    input:
        "/home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{experiment}/{line}"
    output:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/qc/{experiment}_{line}_multiqc_report.html"
    shell:
        "cd {input} && multiqc . -o multiqc && cp multiqc/multiqc_report.html {output}"

rule bwa_map:
    input:
        "/home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta", #ancestor/ancestor_consensus.fasta
        "/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz",
        "/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R2_001.fastq.gz"
    output:
        temp("/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}.bam.temp")
    log:
        "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/logs/{experiment}_{line}_{sample}.log"
    threads: 1
    shell:
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"



rule samtools_sort:
    input:
        "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}.bam.temp"
    output:
        "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam"
    shell:
        "samtools sort {input} > {output}"




rule samtools_index:
    input:
        "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam"
    output:
        "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam.bai"
    shell:
        "samtools index {input}"



rule mapping_quality:
    input:
        bam="/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam"
    output:
        "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted_stats/qualimapReport.html"
    shell:
        "qualimap bamqc -bam {input.bam}"


rule coverage_assessment:
    input:
        bam="/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam",
        bai="/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam.bai"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/coverage-assessment/pipeline-outputs/coverage-plots/{experiment}-{line}-{sample}.png",
        "/home/amovas/data/genome-evo-proj/data/processed-data/coverage-assessment/pipeline-outputs/{experiment}/{line}/{sample}.quals"
    script:
        "/home/amovas/data/genome-evo-proj/src/snakemake-run/python-scripts/quality_control_mapping.py"

rule collect_coverage_assessment:
    input:
        #expand("quality_control/{exp}/{line}/{sample}.quals",
        #       zip, exp=directories,line=lines, sample=IDS),
        expand("/home/amovas/data/genome-evo-proj/data/processed-data/coverage-assessment/pipeline-outputs/{experiment}/{line}/{sample}.quals",
                zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/qc/all_quals.tsv.gz"
    shell:
        "cat {input} | gzip >> {output}"




rule mapping_to_mutations:
    input:
        bam="/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam",
        bai="/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam.bai"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/mutations/pipeline-outputs/{experiment}/{line}/{sample}.csv"
    shell:
        "python /home/amovas/data/genome-evo-proj/src/snakemake-run/python-scripts/my_call_mutations.py {input.bam} /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta NL43_ann_wk0virusPassRef_plasmid > {output}"

# AF324493.2
# NL43seq210314
# K03455.1
# NL43_ann_wk0virusPassRef_SEQregion



rule collect_mutations:
    input:
        expand("/home/amovas/data/genome-evo-proj/data/processed-data/mutations/pipeline-outputs/{experiment}/{line}/{sample}.csv",
                zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        #"mutations/all_ins.pkl",
        #"mutations/all_dels.pkl",
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.csv.gz"
    shell:
        "python /home/amovas/data/genome-evo-proj/src/snakemake-run/python-scripts/my_collect_mutations.py {input}"


rule variant_to_vcf:
    input:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.csv.gz"
    output:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.vcf.gz"
    shell:
        "Rscript /home/amovas/data/genome-evo-proj/src/snakemake-run/R-scripts/variant_to_vcf.R && gzip /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.vcf"



rule annotate_vcf:
    input:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.vcf.gz"
    output:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.ann.vcf.gz"
    shell:
        "cd /home/amovas/data/genome-evo-proj/data/reference/annotations && snpEff -v hiv_plasmid_ref_genome /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.vcf.gz > /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.ann.vcf \
	&& gzip /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.ann.vcf"


rule prep_feature_list:
    input:
        "/home/amovas/data/genome-evo-proj/data/reference/annotations/features/sequence-features.tsv"
    output:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/cds_feature_list.tsv"
    script:
        "/home/amovas/data/genome-evo-proj/src/snakemake-run/R-scripts/prep_feature_list.R"


rule extract_annotation:
    input:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/cds_feature_list.tsv",
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_variants.ann.vcf.gz"
    output:
        "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_annotated_variants.tsv.gz"
    resources:
        mem_mb=20000
    shell:
        "Rscript /home/amovas/data/genome-evo-proj/src/snakemake-run/R-scripts/aa_change.R && Rscript /home/amovas/data/genome-evo-proj/src/snakemake-run/R-scripts/add_mut_context.R"
