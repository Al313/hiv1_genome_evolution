

wd="/home/amovas/data/genome-evo-proj"
EXP,LINES,SAMPLES = glob_wildcards("/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz")



# If you want to run the pipeline for a subset of samples use these assignments:

#EXP = ["iv","iv", "iii", "iii", "iii", "iii"]
#LINES = ["20","20", "15", "15", "15", "15"]
#SAMPLES = ["20MT2EXPIVVP180seq09052019_S24_L001","20MT2EXPIVVP410combinedseq08122022_S33_L001", "15MT4EXPIIIVP30seq08122022_S6_L001", "15MT4EXPIIIVP130seq18042019_S13_L001", "15MT4EXPIIIVP440seq20052022_S16_L001", "15MT4EXPIIIVP510seq09092022_S9_L001"]


rule all:
    input:
       #"/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/all_annotated_variants.tsv.gz",
       "/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/qc/all_quals.tsv.gz",
       #"/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta.iqtree",
       #"/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta.iqtree",
       expand("/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/qc/{experiment}_{line}_multiqc_report.html", zip, experiment=EXP, line=LINES),
       expand("/home/amovas/data/genome-evo-proj/data/processed-data/fastqc-reports/pipeline-outputs/{experiment}/{line}/{sample}_R1_001_fastqc.html", zip,experiment=EXP, line=LINES, sample=SAMPLES),
       expand("/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted_stats/qualimapReport.html", zip,experiment=EXP, line=LINES, sample=SAMPLES),
       #expand("/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_1P.fastq.gz", zip,experiment=EXP, line=LINES, sample=SAMPLES),
       #expand("/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_2P.fastq.gz", zip,experiment=EXP, line=LINES, sample=SAMPLES),
       #expand("/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam", zip,experiment=EXP, line=LINES, sample=SAMPLES)


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
        "cd {input} && multiqc . -o multiqc && cp multiqc/multiqc_report.html {output} && rm -rf multiqc"

rule read_trim:
    input:
        one="/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R1_001.fastq.gz",
        two="/home/amovas/data/genome-evo-proj/data/freezed-raw-data/fastq/{experiment}/{line}/{sample}_R2_001.fastq.gz"
    output:
        one="/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_1P",
        two="/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_2P"
    shell:
        "x={output.one} && y=${{x%%_1P*}} && trimmomatic PE -phred33 {input.one} {input.two} -baseout ${{y}} LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36"

rule read_trim_post_processing:
    input:
        one="/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_1P",
        two="/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_2P"
    output:
        one="/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_1P.fastq.gz",
        two="/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_2P.fastq.gz"
    shell:
        "gzip -S '.fastq.gz' {input.one} && gzip -S '.fastq.gz' {input.two}"

rule bwa_map:
    input:
        "/home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta", #ancestor/ancestor_consensus.fasta
        "/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_1P.fastq.gz",
        "/home/amovas/shared/genome-evo-proj/data/freezed-raw-data/fastq/trimmed/{experiment}/{line}/{sample}_2P.fastq.gz"
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

####

rule get_consensus:
    input:
        "/home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/{experiment}/{line}/{sample}_sorted.bam"
    output:
        tmp1=temp("/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/{experiment}/{line}/{sample}.fa"),
        tmp2=temp("/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/{experiment}/{line}/{sample}.fa1"),
        final_cons="/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/{experiment}/{line}/{sample}.fasta"
    shell:
        "samtools consensus --mode simple -f fasta {input} --show-del yes --show-ins no -aa --call-fract 0.5 -o {output.tmp1} && seqtk seq {output.tmp1} > {output.tmp2} && sed -n 2p {output.tmp2} | cut -c790-2085 > {output.final_cons} && sed -i '1i>{wildcards.sample}' {output.final_cons}"

rule collect_consensus:
    input:
        expand("/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/{experiment}/{line}/{sample}.fasta", zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full.fasta"
    shell:
        "cat {input} >> {output}"


rule consensus_line_swap:
   input:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full.fasta"
   output:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full2.fasta"
   shell:
       "cat {input} | grep -A 1 -E '>(13|14|15|16).*' | sed -E 's/>15(MT4EXPIIIVP[4,5])/>a\\1/g' | sed -E 's/>16(MT4EXPIIIVP[4,5])/>b\\1/g' | sed -E 's/>a(MT4EXPIIIVP[4,5])/>16\\1/g' | sed -E 's/>b(MT4EXPIIIVP[4,5])/>15\\1/g' > {output} && rm {input}"




rule consensus_50th_interval:
   input:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full2.fasta"
   output:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full.fasta"
   shell:
       "cat {input} | grep -A 1 -E '>(13|14|15|16).*VP(50|100|150|200|250|300|350|400|450|500|550|570)[^0-9+\-]' | grep -v -- '^--$' > {output}"


rule consensus_additions_all:
    input:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full2.fasta"
    output:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full3.fasta"
    shell:
       "cat {input} > {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta | cut -c790-2085 >> {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -d ' ' -f1 >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -c790-2085 >> {output}"


rule consensus_additions_50th:
    input:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full.fasta"
    output:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full2.fasta"
    shell:
       "cat {input} > {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta | cut -c790-2085 >> {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -d ' ' -f1 >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -c790-2085 >> {output}"


rule msa_all:
   input:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full3.fasta"
   output:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full_msa.fasta"
   shell:
       "muscle -super5 {input} -output {output}"

rule msa_50th:
   input:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full2.fasta"
   output:
       "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full_msa.fasta"
   shell:
       "muscle -align {input} -output {output}"


rule msa_post_processing_all:
    input:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full_msa.fasta"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta"
    shell:
       "sed -i -e 's/*/-/g' {input} && cat {input} | sed -E 's/([0-9]{{2}})MT[0-9]EXPIIIVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' | sed -E 's/([0-9]{{2}})MT[0-9]EXPIVVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' > {output} && \
        sed -i -e 's/13:/MT-2_1:/g' {output} && sed -i -e 's/14:/MT-2_2:/g' {output} && sed -i -e 's/15:/MT-4_1:/g' {output} && sed -i -e 's/16:/MT-4_2:/g' {output}"



rule msa_post_processing_50th:
    input:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full_msa.fasta"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta"
    shell:
       "sed -i -e 's/*/-/g' {input} && cat {input} | sed -E 's/([0-9]{{2}})MT[0-9]EXPIIIVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' | sed -E 's/([0-9]{{2}})MT[0-9]EXPIVVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' > {output} && \
        sed -i -e 's/13:/MT-2_1:/g' {output} && sed -i -e 's/14:/MT-2_2:/g' {output} && sed -i -e 's/15:/MT-4_1:/g' {output} && sed -i -e 's/16:/MT-4_2:/g' {output}"


rule tree_inference_all:
    input:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta.iqtree"
    shell:
        "iqtree -s {input} -m MFP -B 1000 -redo"


rule tree_inference_50th:
    input:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta"
    output:
        "/home/amovas/data/genome-evo-proj/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta.iqtree"
    shell:
        "iqtree -s {input} -m MFP -B 1000 -redo"


####

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


