
# load libraries

library("stringr")
library("dplyr")
library("magrittr")


# determine the server path
if (file.exists("/home/amovas/")){
  print("Remote HPC Connection!")
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  print("Local PC Connection!")
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
}

# read in raw feature list downloaded as gff3 from NCBI for HIV-1 NL4-3 strain
features <- read.table(file = paste0(wd,"data/reference/annotations/features/sequence-features.tsv"), sep = "\t", stringsAsFactors = F, header = T)

# assign the column names (gff3 format)
colnames(features) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# setting the correct LTR coordination for 5' and 3'
features$end[features$type == "long_terminal_repeat" & features$start == 1] <- 634
features$start[features$type == "long_terminal_repeat" & features$end == 634] <- 551 # from beginning of 5' repeat region till the end of 5' utr
features$end[features$type == "long_terminal_repeat" & features$start == 9076] <- 9528 # from beginning of 3'utr till the end of 3' repeat region





features %<>% mutate(type = ifelse(type=="long_terminal_repeat", "untranslated_region", type))


# filter out features other than CDS
features_of_interest <- features[features$type %in% c("CDS", "untranslated_region"),]


row.names(features_of_interest) <- 1:nrow(features_of_interest)

# assign the ORF 
head(features_of_interest)
features_of_interest$ORF <- NA
features_of_interest$ORF[features_of_interest$type == "CDS"] <- (features_of_interest$start[features_of_interest$type == "CDS"]+as.numeric(features_of_interest$phase[features_of_interest$type == "CDS"])-455+1)%%3

features_of_interest$ORF[features_of_interest$ORF==0]<-3


# extract gene names

attributes <- str_split(features_of_interest$attributes, pattern = ";")
att_of_interest_bool <- sapply(attributes, function(x) str_detect(x,pattern = "product"))
att_of_interest_bool[[1]][2] <- TRUE
att_of_interest_bool[[13]][2] <- TRUE


att_of_interest <- mapply(function(x,y) x[[y]], x=attributes, y=sapply(att_of_interest_bool, "which"))

product_name <- sapply(str_split(sapply(str_split(att_of_interest, pattern = " "), "[[", 1), pattern = "="), "[[",2)

features_of_interest$feature <- product_name
features_of_interest$feature[1] <- "5UTR"
features_of_interest$feature[13] <- "3UTR"

# extract gene ids


attributes <- str_split(features_of_interest$attributes, pattern = ";")
attributes[[1]] <- append(attributes[[13]], "Name=AF324490.5")
attributes[[13]] <- append(attributes[[13]], "Name=AF324499.3")

att_of_interest_bool <- sapply(attributes, function(x) str_detect(x,pattern = "Name"))

att_of_interest <- mapply(function(x,y) x[[y]], x=attributes, y=sapply(att_of_interest_bool, "which"))

feature_id <- sapply(str_split(att_of_interest, pattern = "="), "[[",2)

features_of_interest$feature_id <- feature_id



# add 1-based index range of features to the table

features_of_interest$feature_range <- paste(features_of_interest$feature, features_of_interest$start, features_of_interest$end, sep = "_")

features_of_interest <- features_of_interest[,-c(4,5,8)]

a <- aggregate(feature_range~feature_id,features_of_interest,FUN=paste, collapse = "|")

b <- aggregate(ORF~feature_id,features_of_interest,FUN=paste, collapse = "|", na.action=NULL)

features_of_interest <- merge(merge(features_of_interest[,-c(7,10)], a, by="feature_id"), b, by="feature_id") %>%
    distinct()

# new line characters mess things up in attribute field
features_of_interest <- features_of_interest[,-7]


## add additional features to the list

features_of_interest <- rbind(features_of_interest,c("AF324489.5","AF324493.2","Genbank","untranslated_region",".","+","5R","5R_454_550",NA), c("AF324491.5","AF324493.2","Genbank","untranslated_region",".","+","5LTRLS","5LTRLS_635_789",NA),c("AF324500.3","AF324493.2","Genbank","untranslated_region",".","+","3R","3R_9529_9625",NA)) 


# write the table

write.table(features_of_interest, file = paste0(wd,"results/tables/pipeline-outputs/cds_feature_list.tsv"), sep = "\t", row.names = F, quote = F)

