
# set working directory paths
wd="/home/amovas/data/genome-evo-proj"
wd_shared="/home/amovas/shared/genome-evo-proj"


# obtain wildcrad elements
EXP,LINES,SAMPLES = glob_wildcards(f"{wd}/data/freezed-raw-data/fastq/{{experiment}}/{{line}}/{{sample}}_R1_001.fastq.gz")


# to run the pipeline for a subset of samples use these assignments:
# EXP = ["iv","iv", "iii", "iii", "iii", "iii"]
# LINES = ["20","20", "15", "15", "15", "15"]
# SAMPLES = ["20MT2EXPIVVP180seq09052019_S24_L001","20MT2EXPIVVP410combinedseq08122022_S33_L001", "15MT4EXPIIIVP30seq08122022_S6_L001", "15MT4EXPIIIVP130seq18042019_S13_L001", "15MT4EXPIIIVP440seq20052022_S16_L001", "15MT4EXPIIIVP510seq09092022_S9_L001"]


