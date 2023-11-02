
# load libraries

library(ggplot2)
library(stringr)
library(tidyr)

# read in quality table
quals = read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/processed-data/quality_control/3-p/all_quals.tsv.gz", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(quals) = c("file","average","cover","average_a","average_b","average_c","average_d","average_e","cover_all","mapped", "mapquality")
quals$full_sample_name = as.character(lapply(str_split(lapply(str_split(quals$file, pattern = "/"), "[[", 12), pattern = "_"), "[[", 1))
quals = quals[,-1]


# read in the metadata
meta = read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/metadata/NGS_samples_list_processed_vlast.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
meta = meta[meta$included == TRUE,]



# merge the two dataset
quals_meta <- merge(quals,meta, by = "full_sample_name")


# check merged worked fine
nrow(quals_meta)
quals[!(quals$full_sample_name %in% quals_meta$full_sample_name),]


##### plot

# keep relevant fields
quals <- quals_meta[,c("full_sample_name", "average_a", "average_b", "average_c", "average_d", "average_e", "virus_line_no", "cover_all")]

# reshape the dataset and prepare for plotting
quals_gathered <- gather(quals, key = "amplicon", value = "average_depth", 2:6)

# convert certain fields to factor
quals_gathered$virus_line_no <- factor(quals_gathered$virus_line_no, levels = c(13,14,15,16,17,18,19,20))
quals_gathered$amplicon <- factor(quals_gathered$amplicon, levels = c("average_a", "average_b", "average_c", "average_d", "average_e"))
quals_gathered$full_sample_name <- factor(quals_gathered$full_sample_name)

# check the data type of dataset fields
lapply(quals_gathered,class)

# sort the dataset (optional)
quals_gathered <- quals_gathered[order(quals_gathered$full_sample_name),]

# logorithmize the depth value for better visualization
quals_gathered$average_depth <- log10(quals_gathered$average_depth)


# render the ggplot object

coverage_map <- ggplot(data=quals_gathered, aes(x = amplicon, y = average_depth,group=full_sample_name))+
    geom_line(linetype="dashed", color="blue", size=1.2) +
    geom_point(color="red", size=3) +
    geom_text(data=subset(quals_gathered, cover_all == "False" & average_depth <= 3),
    aes(x=amplicon,y=average_depth,label=full_sample_name),size = 6) +
    labs(title = "Coverage Check!\n", x = "\nAmplicons", y = "Average Depth (log10)\n") +
    scale_x_discrete(labels = c("Amp_A","Amp_B","Amp_C","Amp_D","Amp_E")) +
    scale_y_continuous(limit = c(-1,6), breaks = 0:6) +
    geom_hline(yintercept = 3, linetype="dashed") +
    theme_bw() +
    theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text = element_text(size = 15))



# save the ggplot object
png(filename = "/home/amovas/data/genome-evo-proj/results/figs/png/coverage_map.png",
    width = 1500,
    height = 1800)
coverage_map
dev.off()




# a box plot?






