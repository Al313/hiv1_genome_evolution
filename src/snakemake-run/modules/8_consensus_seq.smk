rule getConsensus:
    input:
        f"{wd_shared}/data/processed-data/mappings/pipeline-outputs/{{experiment}}/{{line}}/{{sample}}_sorted.bam"
    output:
        tmp1=temp(f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fa"),
        tmp2=temp(f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fa1"),
        final_cons=f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fasta"
    shell:
        "samtools consensus --mode simple -f fasta {input} --show-del yes --show-ins no -aa --call-fract 0.5 -o {output.tmp1} && seqtk seq {output.tmp1} > {output.tmp2} && sed -n 2p {output.tmp2} | cut -c454-9625 > {output.final_cons} && sed -i '1i>{wildcards.sample}' {output.final_cons}"

rule collectConsensus:
    input:
        expand(f"{wd}/data/processed-data/consensus/pipeline-outputs/full/{{experiment}}/{{line}}/{{sample}}.fasta", zip, experiment=EXP, line=LINES, sample=SAMPLES)
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full.fasta"
    shell:
        "cat {input} >> {output}"

rule consensusLineSwap:
   input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full.fasta"
   output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full2.fasta"
   shell:
       "cat {input} | grep -A 1 -E '>(13|14|15|16).*' | sed -E 's/>15(MT4EXPIIIVP[4,5,6])/>a\\1/g' | sed -E 's/>16(MT4EXPIIIVP[4,5,6])/>b\\1/g' | sed -E 's/>a(MT4EXPIIIVP[4,5,6])/>16\\1/g' | sed -E 's/>b(MT4EXPIIIVP[4,5,6])/>15\\1/g' | sed -E 's/>15(MT4EXPIIIVP[3][7,8,9])/>a\\1/g' | sed -E 's/>16(MT4EXPIIIVP[3][7,8,9])/>b\\1/g' | sed -E 's/>a(MT4EXPIIIVP[3][7,8,9])/>16\\1/g' | sed -E 's/>b(MT4EXPIIIVP[3][7,8,9])/>15\\1/g' > {output} && rm {input}"


rule consensusPostProcessing:
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

