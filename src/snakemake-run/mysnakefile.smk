

# set working directory paths
wd="/home/amovas/data/genome-evo-proj"
wd_shared="/home/amovas/shared/genome-evo-proj"


# obtain wildcrad elements
EXP,LINES,SAMPLES = glob_wildcards(f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R1_001.fastq.gz")


# to run the pipeline for a subset of samples use these assignments:
#EXP = ["iv","iv", "iii", "iii", "iii", "iii"]
#LINES = ["20","20", "15", "15", "15", "15"]
#SAMPLES = ["20MT2EXPIVVP180seq09052019_S24_L001","20MT2EXPIVVP410combinedseq08122022_S33_L001", "15MT4EXPIIIVP30seq08122022_S6_L001", "15MT4EXPIIIVP130seq18042019_S13_L001", "15MT4EXPIIIVP440seq20052022_S16_L001", "15MT4EXPIIIVP510seq09092022_S9_L001"]


# define the master rule
rule all:
    input:
       #expand(f"{wd}/results/tables/pipeline-outputs/{{experiment}}_variants.csv.gz", experiment=set(EXP))
       #expand(f"{wd}/results/tables/pipeline-outputs/{{experiment}}_variants.vcf.gz", experiment=set(EXP))
       #expand(f"{wd}/results/tables/pipeline-outputs/{{experiment}}_variants.ann.vcf.gz", experiment=set(EXP))
       expand(f"{wd}/results/tables/pipeline-outputs/{{experiment}}_annotated_variants.tsv.gz", experiment=set(EXP))
       #f"{wd}/results/tables/pipeline-outputs/qc/all_quals.tsv.gz",
       #expand(f"{wd}/results/tables/pipeline-outputs/qc/{{experiment}}_{{line}}_multiqc_report.html", zip, experiment=EXP, line=LINES),
       #expand(f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_R1_001_fastqc.html", zip,experiment=EXP, line=LINES, sample=SAMPLES),
       #expand(f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted_stats/qualimapReport.html", zip,experiment=EXP, line=LINES, sample=SAMPLES)


### sequencing quelity

rule sequencing_quality:
    input:
        read1=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R1_001.fastq.gz",
        read2=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R2_001.fastq.gz"
    output:
        f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_R1_001_fastqc.html",
        f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_R2_001_fastqc.html"
    shell:
        "fastqc {input.read1} -o /home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{wildcards.experiment}/{wildcards.line} && fastqc {input.read2} -o /home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{wildcards.experiment}/{wildcards.line}"

rule collect_sequencing_quality:
    input:
        f"{wd}/data/processed-data/fastqc-reports/pipeline-outputs/{{experiment}}/{{line}}"
    output:
        f"{wd}/results/tables/pipeline-outputs/qc/{{experiment}}_{{line}}_multiqc_report.html"
    shell:
        "cd {input} && multiqc . -o multiqc && cp multiqc/multiqc_report.html {output} && rm -rf multiqc"

rule read_trim:
    input:
        one=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R1_001.fastq.gz",
        two=f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R2_001.fastq.gz"
    output:
        one=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P",
        two=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P"
    shell:
        "x={output.one} && y=${{x%%_1P*}} && trimmomatic PE -phred33 {input.one} {input.two} -baseout ${{y}} LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36"

rule read_trim_post_processing:
    input:
        one=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P",
        two=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P"
    output:
        one=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P.fastq.gz",
        two=f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P.fastq.gz"
    shell:
        "gzip -S '.fastq.gz' {input.one} && gzip -S '.fastq.gz' {input.two}"

### !

### mapping

rule bwa_map:
    input:
        f"{wd}/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta", #ancestor/ancestor_consensus.fasta
        f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_1P.fastq.gz",
        f"{wd_shared}/data/freezed-raw-data/fastq/trimmed/{{experiment}}/{{line}}/{{sample}}_2P.fastq.gz"
    output:
        temp(f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.bam.temp")
    log:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/logs/{{experiment}}_{{line}}_{{sample}}.log"
    threads: 1
    shell:
        "(bwa mem -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"



rule samtools_sort:
    input:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.bam.temp"
    output:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    shell:
        "samtools sort {input} > {output}"




rule samtools_index:
    input:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    output:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam.bai"
    shell:
        "samtools index {input}"

### !

### mapping quality

rule mapping_quality:
    input:
        bam=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    output:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted_stats/qualimapReport.html"
    shell:
        "qualimap bamqc -bam {input.bam}"


rule coverage_assessment:
    input:
        bam=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam",
        bai=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam.bai"
    output:
        f"{wd}/data/processed-data/coverage-assessment/pipeline-outputs/coverage-plots/{{experiment}}-{{line}}-{{sample}}.png",
        f"{wd}/data/processed-data/coverage-assessment/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.quals"
    script:
        f"{wd}/src/snakemake-run/python-scripts/quality_control_mapping.py"

rule collect_coverage_assessment:
    input:
        #expand("quality_control/{exp}/{{line}}/{{sample}}.quals",
        #       zip, exp=directories,line=lines, sample=IDS),
        expand(f"{wd}/data/processed-data/coverage-assessment/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.quals",
                zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        f"{wd}/results/tables/pipeline-outputs/qc/all_quals.tsv.gz"
    shell:
        "cat {input} | gzip >> {output}"

### !


### variant calling

rule mapping_to_mutations:
    input:
        bam=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam",
        bai=f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam.bai"
    output:
        f"{wd}/data/processed-data/mutations/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}.csv"
    shell:
        "python /home/amovas/data/genome-evo-proj/src/snakemake-run/python-scripts/my_call_mutations.py {input.bam} /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta NL43_ann_wk0virusPassRef_plasmid > {output}"


# Helper function to get samples for a specific experiment
def get_samples_for_experiment(wildcards):
    exp_samples = []
    for i, exp in enumerate(EXP):
        if exp == wildcards.experiment:
            exp_samples.append(f"{wd}/data/processed-data/mutations/pipeline-outputs/{EXP[i]}/{LINES[i]}/{SAMPLES[i]}.csv")
    return exp_samples


rule collect_mutations:
    input:
        get_samples_for_experiment
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.csv.gz"
    shell:
        "python /home/amovas/data/genome-evo-proj/src/snakemake-run/python-scripts/my_collect_mutations.py {input}"

### !

### vcf generation and annotation

rule variant_to_vcf:
    input:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.csv.gz"
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.vcf.gz"
    shell:
        "Rscript /home/amovas/data/genome-evo-proj/src/snakemake-run/R-scripts/variant_to_vcf.R {wildcards.experiment} && gzip /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.vcf"



rule annotate_vcf:
    input:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.vcf.gz"
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.ann.vcf.gz"
    shell:
        "cd /home/amovas/data/genome-evo-proj/data/reference/annotations && snpEff -v hiv_plasmid_ref_genome /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.vcf.gz > /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.ann.vcf \
	&& gzip /home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/{wildcards.experiment}/{wildcards.experiment}_variants.ann.vcf"

### !

### processing of final variant table

rule prep_feature_list:
    input:
        f"{wd}/data/reference/annotations/features/sequence-features.tsv"
    output:
        f"{wd}/results/tables/pipeline-outputs/cds_feature_list.tsv"
    script:
        f"{wd}/src/snakemake-run/R-scripts/prep_feature_list.R"


rule extract_annotation:
    input:
        f"{wd}/results/tables/pipeline-outputs/cds_feature_list.tsv",
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_variants.ann.vcf.gz"
    output:
        f"{wd}/results/tables/pipeline-outputs/{{experiment}}/{{experiment}}_annotated_variants.tsv.gz"
    resources:
        mem_mb=50000
    shell:
        "Rscript /home/amovas/data/genome-evo-proj/src/snakemake-run/R-scripts/aa_change.R {wildcards.experiment} && Rscript /home/amovas/data/genome-evo-proj/src/snakemake-run/R-scripts/add_mut_context.R {wildcards.experiment}"

### !
