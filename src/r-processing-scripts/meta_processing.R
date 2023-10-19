
# change the name of previous metadata to include version before starting

library(dplyr)
library(stringr)

# set variables to process relevant metadata
arrival_date <- "04082023"
ngs_sample_list_version <- "414"

# read in the relevant metadata
ngs_sample_list_new <- read.csv(file = paste0("/home/amovas/data/genome-evo-proj/data/incoming-raw-data/", arrival_date, "/NGS_samples_list_all_runs_NGS_R1R2_v", ngs_sample_list_version, "_100823.csv"), sep = ",", header = T, stringsAsFactors = F)
ngs_sample_list_new[is.na(ngs_sample_list_new)] <- "NA"



# filter for relevant samples

ngs_sample_list_new <- ngs_sample_list_new[str_detect(ngs_sample_list_new$forward_name, pattern = "VPIII_p520-570"),]

# check if the number of samples match
nrow(ngs_sample_list_new)


# clean the data: retain relevant fields, rename, and re-order to be consistent with my metadata version
relevant_fields <- c("NGS.date", "PROJECT", "NGS.instrument", "NGS.technology", "sample.s.name._.name.for.the.NGS.data", "sample.set.No.", "sample.date", "organism", "sample.source", "amplicon",
"HIV.1.subtype", "forward_name", "reverse_name", "extra_base_uuid", "base_uuid")


ngs_sample_list_new <- ngs_sample_list_new[,relevant_fields]
ngs_sample_list_new$PROJECT <- "HIV_passaging"

colnames(ngs_sample_list_new) <- c("ngs_date", "project", "ngs_instrument", "ngs_technology", "full_sample_name", "sample_set_no", "sample_date", "organism", "sample_source", "amplicon", "organism_subtype", "forward_name", "reverse_name", "extra_base_uuid", "base_uuid")

# extract following info and store them in a new field


ngs_sample_list_new$virus_line_no <- sapply(str_split(ngs_sample_list_new$full_sample_name, "MT"), "[[",1)
ngs_sample_list_new$cell_line_no <- sapply(str_split(ngs_sample_list_new$full_sample_name, ""), "[[",5)

ngs_sample_list_new$exp_no <- sapply(str_split(sapply(str_split(ngs_sample_list_new$full_sample_name, "EXP"), "[[",2), "VP"), "[[", 1)

# make sure number of experiments are consistent
ngs_sample_list_new$exp_no[ngs_sample_list_new$exp_no == "4"] <- "IV"
ngs_sample_list_new$exp_no[ngs_sample_list_new$exp_no == "3"] <- "III"

ngs_sample_list_new$transfer_no <- sapply(str_split(sapply(str_split(ngs_sample_list_new$full_sample_name, "VP"), "[[",2), "seq"), "[[", 1)


# set used,paired_sample,included and update later as needed

ngs_sample_list_new$used <- "T"
ngs_sample_list_new$paired_sample <- NA
ngs_sample_list_new$included <- TRUE


# re-order

ngs_sample_list_new <- ngs_sample_list_new[,c(2, 8, 11, 9, 5, 16:19, 1, 3:4, 6:7, 10, 12:15, 20:22)]

head(ngs_sample_list_new)



# read in the previous metadata
## note that previously to process and clean up the initial meta data I had to run the following code:


# meta$full_sample_name <- str_replace(meta$full_sample_name, pattern = "EXP4", replacement = "EXPIV")
# meta$full_sample_name <- str_replace(meta$full_sample_name, pattern = "EXP3", replacement = "EXPIII")
# meta$full_sample_name <- str_replace(meta$full_sample_name, pattern = "combined", replacement = "")
# meta$full_sample_name[meta$used == "T/2"] <- str_replace(meta$full_sample_name[meta$used == "T/2"], pattern = "seq", replacement = "combinedseq")
# meta$full_sample_name[meta$virus_line_no %in% c(13,14,17:20)] <- str_replace(meta$full_sample_name[meta$virus_line_no %in% c(13,14,17:20)], pattern = "MT4", replacement = "MT2")


ngs_sample_list_old <- read.csv(file = paste0("/home/amovas/data/genome-evo-proj/data/freezed-raw-data/metadata/NGS_samples_list_processed_v0.csv"), sep = ",", header = T, stringsAsFactors = F)


ngs_sample_list_updated <- rbind(ngs_sample_list_old, ngs_sample_list_new)


########### decide how you are gonna handle repeats based on the visual inspection of coverage
########### for 04082023 update

# exclude 16MT4EXPIIIVP510seq09092022 
# combine 13MT2EXPIIIVP450seq20052022 and 13MT2EXPIIIVP450seq13072023


ngs_sample_list_updated[ngs_sample_list_updated$full_sample_name == "16MT4EXPIIIVP510seq09092022",c("used", "included")] <- c("F", FALSE)
ngs_sample_list_updated[ngs_sample_list_updated$full_sample_name == "13MT2EXPIIIVP450seq20052022", c("used", "included")] <- c("T/2", FALSE)
ngs_sample_list_updated[ngs_sample_list_updated$full_sample_name == "13MT2EXPIIIVP450seq13072023", c("full_sample_name", "used", "paired_sample","included")] <- c("13MT2EXPIIIVP450combinedseq13072023", "T/2", "13MT2EXPIIIVP450seq20052022", TRUE)


# prepare the dataframe for saving


ngs_sample_list_updated <- ngs_sample_list_updated %>% mutate_at(c('virus_line_no', 'transfer_no'), as.numeric)


ngs_sample_list_updated <-  ngs_sample_list_updated[with(ngs_sample_list_updated, order(virus_line_no, transfer_no)), ]


# save the updated metadata
write.table(ngs_sample_list_updated, file = "/home/amovas/data/genome-evo-proj/data/freezed-raw-data/metadata/NGS_samples_list_processed_vlast.csv", sep = ",", quote = F, row.names = F)

