
# bash commands:
# 
# for i in *.txt; do paste -d "\t"  - - < ${i} > tab_${i}; done
# rm -v !(tab_*)




# library(tidyr)
# library(stringr)
# 
# 
# ss <- data.frame()
# 
# 
# for (i in 0:5){
#   print(i)
#   
#   for (j in 13:20){
#     print(j)
#     dd <- read.csv(file = paste0("/home/ali313/Desktop/ltee_raw/test-snakemake/freezed-data/covergae/tab_", j, "_", i, ".txt"), header = F, sep = "\t", stringsAsFactor = F)
#     
#     colnames(dd) <- c("transfer_no", paste0("coverage"))
#     dd$region <- i
#     dd$virus_line_no <- j
#     
#     ss <- rbind(ss, dd)
#   }
# }
# 
# # ss[ss$transfer_no == "130"& ss$virus_line_no == 18,]
# 
# 
# ss$transfer_no <- gsub('[^0-9.-]', '', ss$transfer_no)
# 
# ss$virus_line_transfer_no <- paste(ss$virus_line_no, ss$transfer_no, sep = "_")
# 
# ss <- ss[,-c(1,4)]
# 
# ss$region <- paste0("region_", ss$region)
# 
# 
# ss <- spread(ss, region, coverage)
# 
# 
# 
# meta2 <- merge(meta, ss, by = "virus_line_transfer_no", all.x = T)
# 
# # meta$full_sample_name[meta$full_sample_name %notin% meta2$full_sample_name]
# 
# meta2 <- meta2[,-match("virus_line_transfer_no", colnames(meta2))]
# 
# # meta2[rowSums(meta2[,23:28] <= 1000) >= 1,]
# 
# 
# write.table(meta2, file = "/home/ali313/Documents/studies/phd/researcher/data/raw-freezed-data/NGS_samples_list_processed_with_coverage.csv", sep = ",", quote = F, row.names = F)




meta2 <- read.csv(file = "/home/ali313/Documents/studies/phd/researcher/data/raw-freezed-data/NGS_samples_list_processed_with_coverage.csv", sep = ",", header = T, stringsAsFactors = F)




meta2[is.na(meta2)] <- "NA"
meta2[,23:28] <- apply(meta2[,23:28],2, as.numeric)
meta2$whole_virus_coverage <- rowMeans(meta2[,23:28], na.rm = T)
meta2$virus_line_transfer_no <- paste(meta2$virus_line_no, meta2$transfer_no, sep = "_")


# meta2[rowSums(meta2[,23:28] <= 1700) >= 1,]








