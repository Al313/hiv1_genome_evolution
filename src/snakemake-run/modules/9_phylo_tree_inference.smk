rule doMSA:
   input:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full3.fasta"
   output:
       f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full_msa.fasta"
   shell:
       "muscle -super5 {input} -output {output}"

rule MSAPostProcessing:
    input:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/all_consensus_full_msa.fasta"
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta"
    shell:
       "sed -i -e 's/*/-/g' {input} && cat {input} | sed -E 's/([0-9]{{2}})MT[0-9]EXPIIIVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' | sed -E 's/([0-9]{{2}})MT[0-9]EXPIVVP([0-9]{{2,3}}).*_L001?/\\1:\\2/g' > {output} && \
        sed -i -e 's/13:/MT-2_1:/g' {output} && sed -i -e 's/14:/MT-2_2:/g' {output} && sed -i -e 's/15:/MT-4_1:/g' {output} && sed -i -e 's/16:/MT-4_2:/g' {output}"

rule reconstructTree:
    input:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta"
    output:
        f"{wd}/data/processed-data/consensus/pipeline-outputs/full/tree/all_consensus_full_msa2.fasta.iqtree"
    shell:
        "iqtree -s {input} -m MFP -B 1000 -redo"

