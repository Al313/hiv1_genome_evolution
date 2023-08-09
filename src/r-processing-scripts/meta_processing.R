# 
# # data from shared folder (annotated with USED info)
# 
# # meta[meta$patient.s.ID.content == "MT-2_19 EXP4 VP 280 ", "patient.s.ID.content"] <- "MT-2_19 EXP4 VP 280"
# # meta[meta$patient.s.ID.content == "MT-2_20 EXPIII VP 410", "patient.s.ID.content"] <- "MT-2_20 EXPIV VP 410"
# # meta[meta$patient.s.ID.content == "MT-2_20 EXPIII VP 420", "patient.s.ID.content"] <- "MT-2_20 EXPIV VP 420"
# # meta[meta$patient.s.ID.content == "MT-2_19 EXP4 VP 280", "patient.s.ID.content"] <- "MT-2_19 EXPIV VP 280"
# # meta[meta$patient.s.ID.content == "MT-2_17 EXPIII VP 420", "patient.s.ID.content"] <- "MT-2_17 EXPIV VP 420"
# # meta[meta$patient.s.ID.content == "MT-2_18 EXPIII VP 420", "patient.s.ID.content"] <- "MT-2_18 EXPIV VP 420"
# 
# library(dplyr)
# library(stringr)
# 
# 
# meta <- read.csv(file = "/home/ali313/Documents/studies/phd/researcher/data/raw-freezed-data/NGS_samples_list_processed.csv", sep = ",", header = T, stringsAsFactors = F)
# meta <- meta[meta$USED != "",]
# 
# meta[is.na(meta)] <- "NA"
# 
# meta <- meta[,-match(c("patient.s.ID.content", "other.ID", "other.study"), colnames(meta))] #"NGS.run.", "sample.source", "GS.FLX.gasket.s.size", "Order.ID", "ZPHI.ID", "SHCS.ID", "patient.s.initial.or.Lab.ID", "Viral.load.titre", "X", "X.1", "X.2"
# meta$PROJECT <- "HIV_passaging"
# 
# 
# colnames(meta) <- c("ngs_date", "project", "ngs_instrument", "ngs_technology", "full_sample_name", "set", "sample_set_no", "sample_date", "organism", "sample_source", "amplicon", "organism_subtype", "forward_name", "reverse_name", "extra_base_uuid", "base_uuid", "used")
# 
# str_split(meta$full_sample_name[1], "")[[1]][1]
# 
# meta$virus_line_no <- sapply(str_split(meta$full_sample_name, "MT"), "[[",1)
# meta$cell_line_no <- sapply(str_split(meta$full_sample_name, ""), "[[",5)
# 
# meta$exp_no <- sapply(str_split(sapply(str_split(meta$full_sample_name, "EXP"), "[[",2), "VP"), "[[", 1)
# 
# meta$exp_no[meta$exp_no == "4"] <- "IV"
# meta$exp_no[meta$exp_no == "3"] <- "III"
# 
# meta$transfer_no <- sapply(str_split(sapply(str_split(meta$full_sample_name, "VP"), "[[",2), "seq"), "[[", 1)
# meta$transfer_no[str_detect(meta$transfer_no, "Full")] <- sapply(str_split(meta$transfer_no[str_detect(meta$transfer_no, "Full")], "F"), "[[",1)
# 
# 
# 
# 
# meta <- meta[meta$included != "F",]
# 
# 
# 
# 
# head(meta, n = 1)
# 
# 
# 
# 
# 
# # check the merged samples
# unique(paste0(meta$transfer_no[meta$included == "T/2"], meta$virus_line_no[meta$included == "T/2"])) %in% names(table(paste0(meta$transfer_no, meta$virus_line_no))[table(paste0(meta$transfer_no, meta$virus_line_no)) == 2])
# 
# # more chekcs
# 
# # meta[meta$virus_line_no == 13 & meta$cell_line_no == "4",]
# # 
# # 
# # meta$full_sample_name[meta$full_sample_name == "20MT4EXPIVVP10seq22072022"] <- "20MT2EXPIVVP10seq22072022"
# # meta$cell_line_no[meta$full_sample_name == "20MT2EXPIVVP10seq22072022"] <- "2"
# # 
# # 
# # meta$full_sample_name[meta$full_sample_name == "20MT4EXPIVVP380seq22072022"] <- "20MT2EXPIVVP380seq22072022"
# # meta$cell_line_no[meta$full_sample_name == "20MT2EXPIVVP380seq22072022"] <- "2"
# # 
# # meta$full_sample_name[meta$full_sample_name == "13MT4EXPIIIVP440seq08122022"] <- "13MT2EXPIIIVP440seq08122022"
# # meta$cell_line_no[meta$full_sample_name == "13MT2EXPIIIVP440seq08122022"] <- "2"
# 
# 
# 
# 
# for (line in unique(meta$virus_line_no)){
#   print(line)
#   print(seq(10, 510, 10) %in% meta$transfer_no[meta$virus_line_no == line])
# }
# 
# nrow(meta)
# colnames(meta)
# meta <- meta[,c(2, 9, 12, 10, 5, 18:21, 1, 3, 4, 7, 8, 11, 13:16, 17)]
# 
# 
# 
# 
# 
# meta <- meta %>% mutate_at(c('virus_line_no', 'transfer_no'), as.numeric)
# 
# 
# meta <-  meta[with(meta, order(virus_line_no, transfer_no)), ]
# 
# 
# 
# meta$virus_line_transfer_no <- paste(meta$virus_line_no, meta$transfer_no, sep = "_")
# 
# 
# for (i in meta$virus_line_transfer_no[duplicated(meta$virus_line_transfer_no)]){
#   
#   ss <- meta$full_sample_name[meta$virus_line_transfer_no == i]
#   
#   meta$paired_sample[meta$virus_line_transfer_no == i] <- rev(ss)
#   
# }
# meta <- meta[,-match("virus_line_transfer_no", colnames(meta))]
# 
# 
# lastly I add a new column to get rid of the duplicated samples for downstream analyses
# 
# 
# meta$included[meta$used == "T/2" & str_detect(meta$full_sample_name, pattern = "08122022")] <- T
# meta$included[meta$used == "T"] <- T
# 
# meta$included[is.na(meta$included)] <- F
# 
# meta <- meta[meta$included,]
# 
# write.table(meta, file = "/home/ali313/Documents/studies/phd/researcher/data/raw-freezed-data/NGS_samples_list_processed.csv", sep = ",", quote = F, row.names = F)






meta <- read.csv(file = "/home/ali313/Documents/studies/phd/researcher/data/raw-freezed-data/NGS_samples_list_processed.csv", sep = ",", header = T, stringsAsFactors = F)


meta[is.na(meta)] <- "NA"
meta$virus_line_transfer_no <- paste(meta$virus_line_no, meta$transfer_no, sep = "_")

