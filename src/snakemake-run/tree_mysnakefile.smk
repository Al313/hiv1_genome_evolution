

# define the master rule
rule all:
    input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta.iqtree",
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta.iqtree"


####
## gag 790-2084
## full 454-9625
rule get_consensus:
    input:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    output:
        tmp1=temp(f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fa"),
        tmp2=temp(f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fa1"),
        final_cons=f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fasta"
    shell:
        "samtools consensus --mode simple -f fasta {input} --show-del yes --show-ins no -aa --call-fract 0.5 -o {output.tmp1} && seqtk seq {output.tmp1} > {output.tmp2} && sed -n 2p {output.tmp2} | cut -c454-9625 > {output.final_cons} && sed -i '1i>{wildcards.sample}' {output.final_cons}"

rule collect_consensus:
    input:
        expand(f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fasta", zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full.fasta"
    shell:
        "cat {input} >> {output}"


rule consensus_line_swap:
   input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full.fasta"
   output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full2.fasta"
   shell:
       "cat {input} | grep -A 1 -E '>(13|14|15|16).*' | sed -E 's/>15(MT4EXPIIIVP[4,5,6])/>a\\1/g' | sed -E 's/>16(MT4EXPIIIVP[4,5,6])/>b\\1/g' | sed -E 's/>a(MT4EXPIIIVP[4,5,6])/>16\\1/g' | sed -E 's/>b(MT4EXPIIIVP[4,5,6])/>15\\1/g' | sed -E 's/>15(MT4EXPIIIVP[3][7,8,9])/>a\\1/g' | sed -E 's/>16(MT4EXPIIIVP[3][7,8,9])/>b\\1/g' | sed -E 's/>a(MT4EXPIIIVP[3][7,8,9])/>16\\1/g' | sed -E 's/>b(MT4EXPIIIVP[3][7,8,9])/>15\\1/g' > {output} && rm {input}"




rule consensus_50th_interval:
   input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full2.fasta"
   output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full.fasta"
   shell:
       "cat {input} | grep -A 1 -E '>(13|14|15|16).*VP(50|100|150|200|250|300|350|400|450|500|550|570)[^0-9+\-]' | grep -v -- '^--$' > {output}"


rule consensus_additions_all:
    input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full2.fasta"
    output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full3.fasta"
    shell:
       "cat {input} > {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta | cut -c454-9625 >> {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -d ' ' -f1 >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -c454-9625 >> {output}"


rule consensus_additions_50th:
    input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full.fasta"
    output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full2.fasta"
    shell:
       "cat {input} > {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta | cut -c454-9625 >> {output} && \
        sed -n 1p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -d ' ' -f1 >> {output} && \
        sed -n 2p /home/amovas/data/genome-evo-proj/data/reference/hxb2/hxb2_seq_ncbi.fasta | cut -c454-9625 >> {output}"


rule msa_all:
   input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full3.fasta"
   output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full_msa.fasta"
   shell:
       "muscle -super5 {input} -output {output}"

rule msa_50th:
   input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full2.fasta"
   output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full_msa.fasta"
   shell:
       "muscle -align {input} -output {output}"


rule msa_post_processing_all:
    input:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full_msa.fasta"
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta"
    shell:
       "sed -i -e 's/*/-/g' {input} && cat {input} | sed -E 's/([0-9]{{2}})MT[0-9]EXPIIIVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' | sed -E 's/([0-9]{{2}})MT[0-9]EXPIVVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' > {output} && \
        sed -i -e 's/13:/MT-2_1:/g' {output} && sed -i -e 's/14:/MT-2_2:/g' {output} && sed -i -e 's/15:/MT-4_1:/g' {output} && sed -i -e 's/16:/MT-4_2:/g' {output}"



rule msa_post_processing_50th:
    input:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/50th_consensus_full_msa.fasta"
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta"
    shell:
       "sed -i -e 's/*/-/g' {input} && cat {input} | sed -E 's/([0-9]{{2}})MT[0-9]EXPIIIVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' | sed -E 's/([0-9]{{2}})MT[0-9]EXPIVVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' > {output} && \
        sed -i -e 's/13:/MT-2_1:/g' {output} && sed -i -e 's/14:/MT-2_2:/g' {output} && sed -i -e 's/15:/MT-4_1:/g' {output} && sed -i -e 's/16:/MT-4_2:/g' {output}"


rule tree_inference_all:
    input:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta"
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta.iqtree"
    shell:
        "iqtree -s {input} -m MFP -B 1000 -redo"


rule tree_inference_50th:
    input:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta"
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/50th_consensus_full_msa2.fasta.iqtree"
    shell:
        "iqtree -s {input} -m MFP -B 1000 -redo"


####

