
library(tidyr)
library(ggplot2)
library(dplyr)

fraction_threshold <- 0.05

if (fraction_threshold == 0.99){
  type_of_mutation <- "Fixated"
} else if (fraction_threshold == 0.5){
  type_of_mutation <- "Majority"
} else if (fraction_threshold == 0.05) {
  type_of_mutation <- "Minority"
}
 



eva_all_muts <- read.csv(file = "/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/data/HIV-LTE_scripts_for_paper/data/all_muts.csv", sep = ",", header = T, stringsAsFactors = F)
eva_all_muts <- eva_all_muts[eva_all_muts$line %in% c("MT23", "MT24", "MT43", "MT44"),]



eva_all_muts$majority <- ifelse(eva_all_muts$fraction >= fraction_threshold, T, F)
eva_all_muts <- eva_all_muts[eva_all_muts$majority,]


############################################################################
me_all_muts <- read.csv(file = "/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/data/processed-data/mutations/all_muts.csv", sep = ",", header = T, stringsAsFactors = F)
me_all_dels <- read.csv(file = "/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/data/processed-data/mutations/all_dels.csv", sep = ",", header = T, stringsAsFactors = F)

head(me_all_dels)
nrow(me_all_dels)
nrow(me_all_dels[me_all_dels$fraction >= 0.9,])
table(me_all_dels$end[me_all_dels$fraction >= 0.9] - me_all_dels$start[me_all_dels$fraction >= 0.9])

nchar(me_all_dels$delserted[me_all_dels$fraction >= 0.5])


# me_all_muts <- me_all_muts[me_all_muts$line != "MT21" | me_all_muts$passage != "30",]
# me_all_muts <- me_all_muts[me_all_muts$line != "MT22" | me_all_muts$passage != "20",]
# me_all_muts <- me_all_muts[me_all_muts$line != "MT41" | me_all_muts$passage != "30",]
# me_all_muts <- me_all_muts[me_all_muts$line != "MT42" | me_all_muts$passage != "20",]
# 
# me_all_muts <- me_all_muts[me_all_muts$line != "MT21" | me_all_muts$passage != "90",]
# me_all_muts <- me_all_muts[me_all_muts$line != "MT22" | me_all_muts$passage != "90",]
# me_all_muts <- me_all_muts[me_all_muts$line != "MT41" | me_all_muts$passage != "90",]
# me_all_muts <- me_all_muts[me_all_muts$line != "MT42" | me_all_muts$passage != "90",]
# 
# me_all_muts$line <- substr(me_all_muts$line, 1,4)

# write.table(me_all_muts[,-1], file = "/home/ali313/Desktop/ltee_raw/test-snakemake/mutations/all_muts_processed.csv", sep = ",", quote = F, row.names = F)



me_all_muts <- me_all_muts[me_all_muts$passage <= 240,]


me_all_muts$majority <- ifelse(me_all_muts$fraction >= fraction_threshold, T, F)
me_all_muts <- me_all_muts[me_all_muts$majority,]



df_processed <- data.frame()


line_conversion <- c(13, 14, 15, 16)

names(line_conversion) <- c("MT23", "MT24", "MT43", "MT44")

i <- 0
for (line in c("MT23", "MT24", "MT43", "MT44")){
  for (passage in as.character(seq(10, 240, by = 10))){
    i <- i +1
    eva_muts_subset <- nrow(eva_all_muts[eva_all_muts$line == line & eva_all_muts$passage == passage,])
    line_me <- as.vector(line_conversion[line])
    me_muts_subset <- nrow(me_all_muts[me_all_muts$line == line_me & me_all_muts$passage == passage,])
    df_processed <- rbind(df_processed, c(line_me, passage, eva_muts_subset, me_muts_subset))
  }
}


colnames(df_processed) <- c("line", "passage", "Eva's counting", "New counting")



df_processed[,2:4] <- apply(df_processed[,2:4], 2, as.numeric)

tibb_counts <- gather(df_processed, key = "method", value = "counts", -line, -passage)

tibb_counts$line <- factor(tibb_counts$line)
tibb_counts$method <- factor(tibb_counts$method)


p <- ggplot(data = tibb_counts, aes(x = passage, y = counts, color = method)) +
  geom_line() + 
  facet_wrap(~line, scales='free') + 
  ylab(paste0(type_of_mutation, " mutation counts")) +
  ggtitle(paste0("Comparison of number of ", tolower(type_of_mutation), " mutations \n")) +
  scale_x_continuous(breaks = seq(0, 240, by = 40), labels = c("0", "40", "80", "120", "160", "200", "240"), limits = c(0,240)) +
  # scale_y_continuous(breaks = seq(0, 50, by = 10), limits = c(0,50)) +
  labs(caption = paste0(type_of_mutation, " mutations: fraction >= ", fraction_threshold)) +
  theme_bw() +
  theme(plot.title = element_text(size = 35, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 90, size = 15, vjust = 0.6, color = "black"),
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(1,1,1.1,1, "cm")) 


png(filename = paste0("/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/results/figs/", tolower(type_of_mutation), "-comparison-plot.png"), height = 480*2, width = 480*2)
p
dev.off()

tibb_counts[tibb_counts$line == "16" & tibb_counts$passage == 200 & tibb_counts$method == "Eva's counting", ]
tibb_counts[tibb_counts$line == "16" & tibb_counts$passage == 200 & tibb_counts$method == "New counting", ]
table(tibb_counts$method)

summary(tibb_counts)

##########
# plotting coverage



eva_all_muts <- read.csv(file = "/home/ali313/Documents/studies/phd/researcher/data/HIV-LTE_scripts_for_paper/data/all_muts.csv", sep = ",", header = T, stringsAsFactors = F)
eva_all_muts <- eva_all_muts[eva_all_muts$line %in% c("MT23", "MT24", "MT43", "MT44"),]


############################################################################
me_all_muts <- read.csv(file = "/home/ali313/Desktop/ltee_raw/test-snakemake/mutations/all_muts.csv", sep = ",", header = T, stringsAsFactors = F)
me_all_muts <- me_all_muts[me_all_muts$line != "MT21" | me_all_muts$passage != "30",]
me_all_muts <- me_all_muts[me_all_muts$line != "MT22" | me_all_muts$passage != "20",]
me_all_muts <- me_all_muts[me_all_muts$line != "MT41" | me_all_muts$passage != "30",]
me_all_muts <- me_all_muts[me_all_muts$line != "MT42" | me_all_muts$passage != "20",]

me_all_muts <- me_all_muts[me_all_muts$line != "MT21" | me_all_muts$passage != "90",]
me_all_muts <- me_all_muts[me_all_muts$line != "MT22" | me_all_muts$passage != "90",]
me_all_muts <- me_all_muts[me_all_muts$line != "MT41" | me_all_muts$passage != "90",]
me_all_muts <- me_all_muts[me_all_muts$line != "MT42" | me_all_muts$passage != "90",]

me_all_muts$line <- substr(me_all_muts$line, 1,4)


# write.table(me_all_muts[,-1], file = "/home/ali313/Desktop/ltee_raw/test-snakemake/mutations/all_muts_processed.csv", sep = ",", quote = F, row.names = F)



me_all_muts <- me_all_muts[me_all_muts$passage <= 240,]







df_eva <- data.frame()


line_conversio <- c("MT21", "MT22", "MT41", "MT42")

names(line_conversio) <- c("MT23", "MT24", "MT43", "MT44")

i <- 0
for (line in c("MT23", "MT24", "MT43", "MT44")){
  for (passage in as.character(seq(10, 240, by = 10))){
    i <- i +1
    eva_muts_subset <- mean(eva_all_muts$coverage[eva_all_muts$line == line & eva_all_muts$passage == passage])
    line_me <- as.vector(line_conversio[line])
    me_muts_subset <- mean(me_all_muts$coverage[me_all_muts$line == line_me & me_all_muts$passage == passage])
    df_eva <- rbind(df_eva, c(line_me, passage, eva_muts_subset, me_muts_subset))
  }
}


colnames(df_eva) <- c("line", "passage", "Eva's coverage", "New coverage")





df_eva[,2:4] <- apply(df_eva[,2:4], 2, as.numeric)

tibb_coverage <- gather(df_eva, key = "method", value = "coverage", -line, -passage)

tibb_coverage$line <- factor(tibb_coverage$line)
tibb_coverage$method <- factor(tibb_coverage$method)

tibb_coverage[tibb_coverage$method == "Eva's coverage" & tibb_coverage$line == "MT42",]
tibb_coverage[tibb_coverage$method == "New coverage" & tibb_coverage$line == "MT42",]


tibb_coverage$coverage[is.nan(tibb_coverage$coverage)] <- 0

tibb_coverage[tibb_coverage$method == "Eva's coverage" & tibb_coverage$line == "MT22" & tibb_coverage$passage == 200,]
tibb_coverage[tibb_coverage$method == "New coverage" & tibb_coverage$line == "MT22" & tibb_coverage$passage == 200,]


agg_tbl <- tibb_coverage %>% group_by(line) %>% 
  summarise(mean_coverage=mean(coverage, na.rm = T),
            .groups = 'drop')


pp <- ggplot(data = tibb_coverage, aes(x = passage, y = coverage, color = method)) +
  geom_line() + 
  facet_wrap(. ~ line, scales = "free") + 
  ylab("sequencing coverage") +
  ggtitle("Comparison of coverage of sequencing \n") +
  scale_x_continuous(breaks = seq(0, 240, by = 40), labels = c("0", "40", "80", "120", "160", "200", "240"), limits = c(0,240)) +
  scale_y_continuous(breaks = seq(0, 40000, by = 10000), labels = c("0", "10k", "20k", "30k", "40k"), limits = c(0,40000)) +
  geom_hline(data = agg_tbl, aes(yintercept = mean_coverage), linetype="dashed", color = "gray") +
  theme_bw() + 
  theme(plot.title = element_text(size = 35, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 90, size = 15, vjust = 0.6, color = "black"),
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(1,1,1.1,1, "cm")) 



png(filename = paste0("/home/ali313/Desktop/ltee_raw/test-snakemake/results/figures/comparison/coverage-comparison-plot.png"), height = 480*2, width = 480*2)
pp
dev.off()










########################################
# processing new data up until the passage #510



# load packages



me_all_muts <- read.csv(file = "/home/ali313/Desktop/ltee_raw/test-snakemake/freezed-data/mutations/all_muts.csv", sep = ",", header = T, stringsAsFactors = F)

exp <- "EXPIII"





me_all_muts_iii <- me_all_muts[me_all_muts$line %notin% c("17", "18", "19", "20"),]

# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "13" | me_all_muts_iii$passage != "30",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "14" | me_all_muts_iii$passage != "20",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "15" | me_all_muts_iii$passage != "30",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "16" | me_all_muts_iii$passage != "20",]
# 
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "13" | me_all_muts_iii$passage != "90",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "14" | me_all_muts_iii$passage != "90",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "15" | me_all_muts_iii$passage != "90",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "16" | me_all_muts_iii$passage != "90",]
# 
# 
# 
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "16" | me_all_muts_iii$passage != "400",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "16" | me_all_muts_iii$passage != "410",]
# 
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "13" | me_all_muts_iii$passage != "420",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "15" | me_all_muts_iii$passage != "420",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "16" | me_all_muts_iii$passage != "420",]
# 
# 
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "14" | me_all_muts_iii$passage != "430",]
# 
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "13" | me_all_muts_iii$passage != "440",]
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "14" | me_all_muts_iii$passage != "440",]
# 
# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$line != "16" | me_all_muts_iii$passage != "480",]
# 
# 
# 
# me_all_muts_iii$line <- substr(me_all_muts_iii$line, 1,2)

me_all_muts_iii$line[me_all_muts_iii$line == 15 & me_all_muts_iii$passage > 360] <- "A"
me_all_muts_iii$line[me_all_muts_iii$line == 16 & me_all_muts_iii$passage > 360] <- "B"

me_all_muts_iii$line[me_all_muts_iii$line == "A"] <- "16"
me_all_muts_iii$line[me_all_muts_iii$line == "B"] <- "15"



# me_all_muts_iii <- me_all_muts_iii[me_all_muts_iii$fraction> 0.05,]










fraction_threshold <- 0.05

if (fraction_threshold == 0.99){
  type_of_mutation <- "Fixated"
} else if (fraction_threshold == 0.5){
  type_of_mutation <- "Majority"
} else if (fraction_threshold == 0.05) {
  type_of_mutation <- "Minority"
}


if (type_of_mutation == "Minority"){
  me_all_muts_iii$majority <- ifelse(me_all_muts_iii$fraction >= fraction_threshold & me_all_muts_iii$fraction < 0.5, T, F)
  print(type_of_mutation)
} else {
  me_all_muts_iii$majority <- ifelse(me_all_muts_iii$fraction >= fraction_threshold, T, F)
  
}

me_all_muts_iii_subset <- me_all_muts_iii[me_all_muts_iii$majority,]






df_processed <- data.frame()


# line_conversio <- c("MT2_1", "MT2_2", "MT4_1", "MT4_2")
# names(line_conversio) <- c("17", "18", "19", "20")
# 


i <- 0
df_processed <- data.frame()
for (line in c("13", "14", "15", "16")){
  for (passage in as.character(seq(10, 510, by = 10))){
    i <- i +1
    line_me <- line
    me_muts_subset <- nrow(me_all_muts_iii_subset[me_all_muts_iii_subset$line == line_me & me_all_muts_iii_subset$passage == passage,])
    df_processed <- rbind(df_processed, c(line_me, passage, me_muts_subset))
  }
}


colnames(df_processed) <- c("line", "passage", "counting")


df_processed[,2:3] <- apply(df_processed[,2:3], 2, as.numeric)

df_processed$line <- factor(df_processed$line)



q <- ggplot(data = df_processed, aes(x = passage, y = counting, color = line)) +
  geom_line(size = 2) + 
  # facet_wrap(~line, scales='free') + 
  ylab(paste0(type_of_mutation, " mutation counts \n")) +
  xlab("\n Passage") +
  ggtitle(paste0("Number of ", tolower(type_of_mutation), " mutations in EXPIII \n")) +
  scale_x_continuous(breaks = seq(0, 510, by = 60), labels = c("0", "60", "120", "180", "240", "300", "360", "420", "480"), limits = c(0,510)) +
  scale_color_discrete(name = "Line", labels=c('MT2_1', 'MT2_2', 'MT4_1', 'MT4_2')) +
  labs(caption = paste0(type_of_mutation, " mutations: fraction >= ", fraction_threshold)) +
  theme_bw() +
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 90, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        plot.caption = element_text(size = 15),
        legend.title=element_text(size=30),
        legend.text = element_text(size=25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(1,1,1.1,1, "cm"))


png(filename = paste0("/home/ali313/Desktop/ltee_raw/test-snakemake/results/figures/new/", exp, "-", tolower(type_of_mutation), "-plot.png"), height = 480*2, width = 480*2)
q
dev.off()





##########
# plotting coverage



exp <- "EXPIII"



df_coverage <- data.frame()


# line_conversio <- c("MT21", "MT22", "MT41", "MT42")
# names(line_conversio) <- c("MT23", "MT24", "MT43", "MT44")

i <- 0
for (line in c("13", "14", "15", "16")){
  for (passage in as.character(seq(10, 510, by = 10))){
    i <- i +1
    line_me <- line
    me_muts_subset <- mean(me_all_muts_iii$coverage[me_all_muts_iii$line == line_me & me_all_muts_iii$passage == passage])
    df_coverage <- rbind(df_coverage, c(line_me, passage, me_muts_subset))
  }
}


colnames(df_coverage) <- c("line", "passage", "coverage")





df_coverage[,2:3] <- apply(df_coverage[,2:3], 2, as.numeric)


df_coverage$line <- factor(df_coverage$line)

df_coverage$coverage[is.nan(df_coverage$coverage)] <- 0

agg_tbl <- df_coverage %>% group_by(line) %>% 
  summarise(mean_coverage=mean(coverage, na.rm = T),
            .groups = 'drop')
agg_df <- as.data.frame(agg_tbl)



qq <- ggplot(data = df_coverage, aes(x = passage, y = coverage, color = line)) +
  geom_line(size = 1) + 
  xlab("Passage") +
  ylab("Average sequencing coverage \n") +
  ggtitle(paste0("Coverage of sequencing data of ", exp, " \n")) +
  scale_x_continuous(breaks = seq(0, 510, by = 60), labels = c("0", "60", "120", "180", "240", "300", "360", "420", "480"), limits = c(0,510)) +
  scale_y_continuous(breaks = seq(0, 60000, by = 10000), labels = c("0", "10k", "20k", "30k", "40k", "50k", "60k"), limits = c(0,60000)) +
  scale_color_discrete(name = "Line", labels=c('MT2_1', 'MT2_2', 'MT4_1', 'MT4_2')) +
  geom_hline(data = agg_df, aes(yintercept = mean_coverage, color = line), linetype="dashed") +
  theme_bw() + 
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 90, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        # plot.caption = element_text(size = 15),
        legend.title=element_text(size=30),
        legend.text = element_text(size=25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(1,1,1.1,1, "cm"))

png(filename = paste0("/home/ali313/Desktop/ltee_raw/test-snakemake/results/figures/new/", exp, "-coverage-plot.png"), height = 480*2, width = 480*2)
qq
dev.off()

