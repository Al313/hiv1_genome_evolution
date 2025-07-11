

# load the libraries

library(Biostrings)
library(stringr)

# get cl arguments
args <- commandArgs(trailingOnly = TRUE)

# Access the first argument
exp <- args[1]


variants_ann <- read.table(file = paste0("/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/", exp, "_annotated_variants.tsv.gz"), 
sep = "\t", stringsAsFactors = FALSE, header = TRUE)


# read in fasta file to retrieve ref information for adjacent positions
fastaFile <- readDNAStringSet("/home/amovas/data/genome-evo-proj/data/reference/plasmid/plasmid-consensus/hiv_plasmid_consensus_genome.fasta")
seq_name = names(fastaFile)
sequence = paste(fastaFile)


# define functions that will extract the mutation contexts
get_mut_context <- function(df, seq){
    pos <- as.numeric(df[2])
    return(str_sub(seq, pos-1,pos+1))
}


get_mut_3_context <- function(df, seq){
    pos <- as.numeric(df[2])
    return(str_sub(seq, pos-2,pos))
}

get_mut_5_context <- function(df, seq){
    pos <- as.numeric(df[2])
    return(str_sub(seq, pos,pos+2))
}

# extract mutation contexts
variants_ann$mut_context <- apply(variants_ann,1,get_mut_context, sequence)
variants_ann$mut_3_context <- apply(variants_ann,1,get_mut_3_context, sequence)
variants_ann$mut_5_context <- apply(variants_ann,1,get_mut_5_context, sequence)


# save results
write.table(variants_ann, 
    file = gzfile(paste0("/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/", exp, "_annotated_variants.tsv.gz")), sep = "\t", row.names = F, quote = F)
