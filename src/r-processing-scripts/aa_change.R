
# load libraries

library(stringr)

# read in raw feature list doownloaded as gff3 from NCBI for HIV-1 NL4-3 strain
features <- read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/reference/annotations/sequence-features.tsv", sep = "\t", stringsAsFactors = F, header = T)

# assign the column names (gff3 format)
colnames(features) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# filter out features other than CDS
features_cds <- features[features$type == "CDS",]
row.names(features_cds) <- 1:nrow(features_cds)

# assign the ORF 
features_cds$ORF <- features_cds$start%%3

# extract gene names

attributes <- str_split(features_cds$attributes, pattern = ";")
att_of_interest_bool <- sapply(attributes, function(x) str_detect(x,pattern = "product"))

att_of_interest <- mapply(function(x,y) x[[y]], x=attributes, y=sapply(att_of_interest_bool, "which"))

product_name <- sapply(str_split(sapply(str_split(att_of_interest, pattern = " "), "[[", 1), pattern = "="), "[[",2)

features_cds$gene <- product_name

View(features_cds)



find.package("stats")




# read annotated variants which is in vcf format

library(vcfR)

ann_var_vcf <- read.vcfR(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/tables/3-p/all_variants.ann.vcf")

ann_var_info <- vcfR::INFO2df(ann_var_vcf)


head(ann_var_info, n = 1)

ann_var_df <- read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/tables/3-p/all_variants.ann.vcf", header = F,
                          sep = "\t", stringsAsFactors = F)


colnames(ann_var_df) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")




