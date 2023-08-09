


# load packages




library(ggplot2)
library(tibble)
library(dplyr)

me_all_muts <- read.csv(file = "/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/data/processed-data/mutations/all_muts.csv", sep = ",", header = T, stringsAsFactors = F)

me_all_muts_iv <- me_all_muts[me_all_muts$line %in% c("17", "18", "19", "20"),]
exp <- "EXPIV"


me_all_muts_iv$line[me_all_muts_iv$line == 18 & me_all_muts_iv$passage > 350] <- "A"
me_all_muts_iv$line[me_all_muts_iv$line == 19 & me_all_muts_iv$passage > 350] <- "B"

me_all_muts_iv$line[me_all_muts_iv$line == "A"] <- "19"
me_all_muts_iv$line[me_all_muts_iv$line == "B"] <- "18"

fraction_thresholds <- c(0.05, 0.5, 0.99)

for (fraction_threshold in fraction_thresholds){
  
  # fraction_threshold <- 0.99
  if (fraction_threshold == 0.99){
    type_of_mutation <- "Fixated"
  } else if (fraction_threshold == 0.5){
    type_of_mutation <- "Majority"
  } else if (fraction_threshold == 0.05) {
    type_of_mutation <- "Minority"
  }
  
  
  if (type_of_mutation == "Minority"){
    me_all_muts_iv$majority <- ifelse(me_all_muts_iv$fraction >= fraction_threshold & me_all_muts_iv$fraction < 0.5, T, F)
  } else {
    me_all_muts_iv$majority <- ifelse(me_all_muts_iv$fraction >= fraction_threshold, T, F)
    
  }
  
  me_all_muts_iv_subset <- me_all_muts_iv[me_all_muts_iv$majority,]
  
  
  
  
  
  
  df_processed <- data.frame()
  
  
  
  
  i <- 0
  df_processed <- data.frame()
  for (line in c("17", "18", "19", "20")){
    for (passage in as.character(seq(10, 510, by = 10))){
      i <- i +1
      line_me <- line
      me_muts_subset <- nrow(me_all_muts_iv_subset[me_all_muts_iv_subset$line == line_me & me_all_muts_iv_subset$passage == passage,])
      df_processed <- rbind(df_processed, c(line_me, passage, me_muts_subset))
    }
  }
  
  
  colnames(df_processed) <- c("line", "passage", "counting")
  
  
  df_processed[,2:3] <- apply(df_processed[,2:3], 2, as.numeric)
  
  
  
  df_processed$line <- factor(df_processed$line)
  
  
  q <- ggplot(data = df_processed, aes(x = passage, y = counting, color = line)) +
    geom_line(size = 2) + 
    ylab(paste0(type_of_mutation, " mutation counts \n")) +
    xlab("\n Passage") +
    ggtitle(paste0("Number of ", tolower(type_of_mutation), " mutations in EXPIV \n")) +
    scale_x_continuous(breaks = seq(0, 510, by = 60), labels = c("0", "60", "120", "180", "240", "300", "360", "420", "480"), limits = c(0,510)) +
    scale_color_discrete(name = "Line", labels=c('MT2_i', 'MT2_ii', 'MT2_iii', 'MT2_iv')) +
    # labs(caption = paste0(type_of_mutation, " mutations: fraction >= ", fraction_threshold)) +
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
  
  
  png(filename = paste0("/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/results/figs/", exp, "-", tolower(type_of_mutation), "-plot.png"), height = 480*2, width = 480*2)
  print(q)
  dev.off()

}

##########
# plotting coverage


exp <- "EXPIV"


df_coverage <- data.frame()


# line_conversio <- c("MT21", "MT22", "MT41", "MT42")
# names(line_conversio) <- c("MT23", "MT24", "MT43", "MT44")

i <- 0
for (line in c("17", "18", "19", "20")){
  for (passage in as.character(seq(10, 510, by = 10))){
    i <- i +1
    line_me <- line
    me_muts_subset <- mean(me_all_muts_iv$coverage[me_all_muts_iv$line == line_me & me_all_muts_iv$passage == passage])
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

# df_coverage[df_coverage$line == 'MT41' & df_coverage$passage == 320,]
# df_coverage[df_coverage$line == 'MT42' & df_coverage$passage == 320,]
# df_coverage[df_coverage$line == 'MT21' & df_coverage$passage == 320,]
# df_coverage[df_coverage$line == 'MT22' & df_coverage$passage == 320,]




qq <- ggplot(data = df_coverage, aes(x = passage, y = coverage, color = line)) +
  geom_line(size = 1) + 
  xlab("Passage") +
  ylab("Average sequencing coverage \n") +
  ggtitle(paste0("Coverage of sequencing data of ", exp, " \n")) +
  scale_x_continuous(breaks = seq(0, 510, by = 60), labels = c("0", "60", "120", "180", "240", "300", "360", "420", "480"), limits = c(0,510)) +
  scale_y_continuous(breaks = seq(0, 60000, by = 10000), labels = c("0", "10k", "20k", "30k", "40k", "50k", "60k"), limits = c(0,60000)) +
  scale_color_discrete(name = "Line", labels=c('MT2_i', 'MT2_ii', 'MT2_iii', 'MT2_iv')) +
  geom_hline(data = agg_df, aes(yintercept = mean_coverage, color = line), linetype="dashed") +
  theme_bw() + 
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 90, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        legend.title=element_text(size=30),
        legend.text = element_text(size=25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(1,1,1.1,1, "cm"))



png(filename = paste0("/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/results/figs/", exp, "-coverage-plot.png"), height = 480*2, width = 480*2)
qq
dev.off()





##########
# clonality (a.k.a intra-population diversity)

me_all_muts <- read.csv(file = "/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/data/processed-data/mutations/all_muts.csv", sep = ",", header = T, stringsAsFactors = F)


me_all_muts_iv <- me_all_muts[me_all_muts$line %in% c("17", "18", "19", "20"),]


me_all_muts_iv$line[me_all_muts_iv$line == 18 & me_all_muts_iv$passage > 350] <- "A"
me_all_muts_iv$line[me_all_muts_iv$line == 19 & me_all_muts_iv$passage > 350] <- "B"

me_all_muts_iv$line[me_all_muts_iv$line == "A"] <- "19"
me_all_muts_iv$line[me_all_muts_iv$line == "B"] <- "18"


me_all_muts_iv_d <- me_all_muts_iv[me_all_muts_iv$fraction> 0.05,]

# d <- density(me_all_muts_iv$fraction)
# plot(d)


d <- ggplot(me_all_muts_iv_d, aes(x = fraction))+
  geom_density(size = 1) +
  ylab(paste0("Density \n")) +
  xlab("Fraction") +
  ggtitle(paste0("Mutation fraction distribution EXPIV \n")) +
  theme_bw() +
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 90, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(1,1,1.1,1, "cm"))


png(filename = paste0("/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/results/figs/EXPIV-fraction-plot.png"), height = 480*2, width = 480*2)
d
dev.off()




# 0.99
# 0.8
# 0.5

exp <- "EXPIV"
min_fraction_threshold <- 0.05


fraction_thresholds2 <- c(0.5, 0.8, 0.99)

for (fraction_threshold in fraction_thresholds2){
  
  # fraction_threshold <- 0.99
  type_of_mutation <- fraction_threshold
  
  
  
  
  i <- 0
  df_processed <- data.frame()
  for (line in c("17", "18", "19", "20")){
    for (passage in as.character(seq(10, 510, by = 10))){
      i <- i +1
      line_me <- line
      clonality_fraction <- nrow(me_all_muts_iv[me_all_muts_iv$line == line_me & me_all_muts_iv$passage == passage & me_all_muts_iv$fraction >= fraction_threshold,])/nrow(me_all_muts_iv[me_all_muts_iv$line == line_me & me_all_muts_iv$passage == passage & me_all_muts_iv$fraction >= min_fraction_threshold,])
      df_processed <- rbind(df_processed, c(line_me, passage, clonality_fraction))
    }
  }
  
  
  colnames(df_processed) <- c("line", "passage", "clonality")
  
  
  df_processed[,2:3] <- apply(df_processed[,2:3], 2, as.numeric)
  
  df_processed$line <- factor(df_processed$line)
  
  df_processed$clonality[is.nan(df_processed$clonality)] <- 0
  
  
  # if (T){
  #   df_processed$clonality[df_processed$line == "18" & df_processed$passage == 130] <- 0.9
  #   df_processed$clonality[df_processed$line == "20" & df_processed$passage == 400] <- 0.9
  # }
  
  
  
  
  q <- ggplot(data = df_processed, aes(x = passage, y = clonality, color = line)) +
    geom_line(size = 2) + 
    xlab("\n Passage") +
    ylab(paste0("Clonality \n")) +
    ggtitle(paste0("HIV-1 population homogeneity in EXPIV \n")) +
    scale_x_continuous(breaks = seq(0, 510, by = 60), labels = c("0", "60", "120", "180", "240", "300", "360", "420", "480"), limits = c(0,510)) +
    scale_color_discrete(name = "Line", labels=c('MT2_i', 'MT2_ii', 'MT2_iii', 'MT2_iv')) +
    # labs(caption = paste0("Number of mutations with fraction >= ", fraction_threshold, "\n _________________________________", "\n Number of mutations with fraction >=", min_fraction_threshold)) +
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
  
  
  
  png(filename = paste0("/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/results/figs/", exp, "-", fraction_threshold, "-diversity-plot.png"), height = 480*2, width = 480*2)
  print(q)
  dev.off()

}


# ==============================================================================

### Exp III

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


# d <- density(me_all_muts_iii$fraction[me_all_muts_iii$fraction> 0.05])
# plot(d, main="Mutation fraction distribution EXPIII")
# polygon(d, col="red", border="blue")


me_all_muts_iii_d <- me_all_muts_iii[me_all_muts_iii$fraction> 0.05,]

# d <- density(me_all_muts_iv$fraction)
# plot(d)


d <- ggplot(me_all_muts_iii_d, aes(x = fraction))+
  geom_density(size = 1) +
  ylab(paste0("Density \n")) +
  xlab("Fraction") +
  ggtitle(paste0("Mutation fraction distribution EXPIII \n")) +
  theme_bw() +
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 90, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(1,1,1.1,1, "cm"))

png(filename = paste0("/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/results/figs/EXPIII-fraction-plot.png"), height = 480*2, width = 480*2)
d
dev.off()


# ------------------------------------------------------------------------------


min_fraction_threshold <- 0.05

fraction_thresholds2 <- c(0.5, 0.8, 0.99)

for (fraction_threshold in fraction_thresholds2){
  
  
  type_of_mutation <- fraction_threshold
  
  i <- 0
  df_processed <- data.frame()
  for (line in c("13", "14", "15", "16")){
    for (passage in as.character(seq(10, 510, by = 10))){
      i <- i +1
      line_me <- line
      clonality_fraction <- nrow(me_all_muts_iii[me_all_muts_iii$line == line_me & me_all_muts_iii$passage == passage & me_all_muts_iii$fraction >= fraction_threshold,])/nrow(me_all_muts_iii[me_all_muts_iii$line == line_me & me_all_muts_iii$passage == passage & me_all_muts_iii$fraction >= min_fraction_threshold,])
      df_processed <- rbind(df_processed, c(line_me, passage, clonality_fraction))
    }
  }
  
  
  colnames(df_processed) <- c("line", "passage", "clonality")
  
  
  df_processed[,2:3] <- apply(df_processed[,2:3], 2, as.numeric)
  
  df_processed$line <- factor(df_processed$line)
  
  df_processed$clonality[is.nan(df_processed$clonality)] <- 0
  
  
  
  
  q <- ggplot(data = df_processed, aes(x = passage, y = clonality, color = line)) +
    geom_line(size = 2) + 
    xlab("\n Passage") +
    ylab(paste0("Clonality \n")) +
    ggtitle(paste0("HIV-1 population homogeneity in EXPIII \n")) +
    scale_x_continuous(breaks = seq(0, 510, by = 60), labels = c("0", "60", "120", "180", "240", "300", "360", "420", "480"), limits = c(0,510)) +
    scale_color_discrete(name = "Line", labels=c('MT2_1', 'MT2_2', 'MT4_1', 'MT4_2')) +
    # labs(caption = paste0("Number of mutations with fraction >= ", fraction_threshold, "\n _________________________________", "\n Number of mutations with fraction >=", min_fraction_threshold)) +
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
  
  
  png(filename = paste0("/Users/alimos313/Documents/studies/phd/research/seq-genome-proj/results/figs/", exp, "-", fraction_threshold, "-diversity-plot.png"), height = 480*2, width = 480*2)
  print(q)
  dev.off()

}


