
# load libraries

library(stringr)

# read in raw feature list doownloaded as gff3 from NCBI for HIV-1 NL4-3 strain
features <- read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/reference/annotations/features/sequence-features.tsv", sep = "\t", stringsAsFactors = F, header = T)


# assign the column names (gff3 format)
colnames(features) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# filter out features other than CDS
features_cds <- features[features$type == "CDS",]
row.names(features_cds) <- 1:nrow(features_cds)

# assign the ORF 

features_cds$ORF <- (features_cds$start+as.numeric(features_cds$phase)-455+1)%%3

features_cds$ORF[features_cds$ORF==0]<-3


# extract gene names

attributes <- str_split(features_cds$attributes, pattern = ";")
att_of_interest_bool <- sapply(attributes, function(x) str_detect(x,pattern = "product"))

att_of_interest <- mapply(function(x,y) x[[y]], x=attributes, y=sapply(att_of_interest_bool, "which"))

product_name <- sapply(str_split(sapply(str_split(att_of_interest, pattern = " "), "[[", 1), pattern = "="), "[[",2)

features_cds$gene <- product_name


# extract gene ids


attributes <- str_split(features_cds$attributes, pattern = ";")
att_of_interest_bool <- sapply(attributes, function(x) str_detect(x,pattern = "Name"))

att_of_interest <- mapply(function(x,y) x[[y]], x=attributes, y=sapply(att_of_interest_bool, "which"))

gene_id <- sapply(str_split(att_of_interest, pattern = "="), "[[",2)

features_cds$gene_id <- gene_id



# add 1-based index range of genes to the table

features_cds$gene_range <- paste(features_cds$gene, features_cds$start, features_cds$end, sep = "_")

features_cds <- features_cds[,-c(4,5,8)]

a <- aggregate(gene_range~gene_id,features_cds,FUN=paste, collapse = "|")

b <- aggregate(ORF~gene_id,features_cds,FUN=paste, collapse = "|")

features_cds <- merge(merge(features_cds[,-c(7,10)], a, by="gene_id"), b, by="gene_id") %>%
    distinct()



# write the table

write.table(features_cds, file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/tables/3-p/cds-feature-list.tsv", sep = "\t", row.names = F, quote = F)

