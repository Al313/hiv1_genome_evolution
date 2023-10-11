

# load the libraries

library(Biostrings)
library(stringr)

variants_ann <- read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/tables/3-p/all_annotated_variants.tsv", 
sep = "\t", stringsAsFactors = FALSE, header = TRUE)


# read in fasta file to retrieve ref information for adjacent positions
fastaFile <- readDNAStringSet("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/reference/plasmid/hiv_plasmid_ref_genome.fasta")
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
    file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/tables/3-p/all_annotated_variants.tsv", sep = "\t", row.names = F, quote = F)
