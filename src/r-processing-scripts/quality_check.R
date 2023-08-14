
library(ggplot2)
library(stringr)

quals = read.table(file = "/home/amovas/data/genome-evo-proj/data/processed-data/quality_control/3-p/all_quals.tsv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(quals) = c("file","average","cover","average_a","average_b","average_c","average_d","average_e","cover_all","mapped", "mapquality")
quals$full_sample_name = as.character(lapply(str_split(lapply(str_split(quals$file, pattern = "/"), "[[", 12), pattern = "_"), "[[", 1))
quals = quals[,-1]
head(quals)


meta = read.table(file = "/home/amovas/data/genome-evo-proj/data/freezed-raw-data/metadata/NGS_samples_list_processed_vlast.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
meta = meta[meta$included == TRUE,]




quals_meta <- merge(quals,meta, by = "full_sample_name")


nrow(quals_meta)
quals[!(quals$full_sample_name %in% quals_meta$full_sample_name),]




