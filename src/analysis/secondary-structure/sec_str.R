

# load libraries

library(readxl)
library(dplyr)

# read-in SHAPE data

shape <- as.data.frame(read_excel("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/external/secondary-structure/41592_2014_BFnmeth3029_MOESM232_ESM.xlsx"))

shape$`1M7 SHAPE MaP`[shape$`1M7 SHAPE MaP` == -999] <- NA


shape$med_55_nt <- NA
for (i in 1:nrow(shape)){
    shape$med_55_nt[i] <- median(shape$`1M7 SHAPE MaP`[max(1,i-27):min(nrow(shape),i+27)], na.rm = T) - median(shape$`1M7 SHAPE MaP`, na.rm = T)
}


plt <- shape %>% ggplot(aes(x = Nt, y = med_55_nt)) +
    geom_line()


pdf("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/shape_reactivity.pdf", width = 21, height = 7)
plt
dev.off()



# read-in base-pairing probabilities

bp_pr <- read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/external/secondary-structure/41592_2014_BFnmeth3029_MOESM234_ESM.txt", sep = "\t", header = T, skip = 1)

bp_pr$pr <- 10**(-bp_pr$X.log10.Probability.)


for (i in unique(bp_pr$i)){
    df <- bp_pr[bp_pr$i == i,]
    
    bp_pr$shannon[bp_pr$i == i] <- -sum(df$pr*(log10(df$pr)))
}





df <- distinct(bp_pr[,c(1,5)])
nrow(df)
View(df)
df$shannon

for (i in 1:9173){
    if (i %notin% df$i){
        df <- rbind(df,c(i,median(df$shannon[df$i %in% max(1,i-5):min(9173,i+5)])))
    }
}






dff <- data.frame()
for (i in 1:9173){
    dff <- rbind(dff, c(i, median(df$shannon[df$i %in% max(1,i-55):min(9173,i+55)], na.rm = T)))
    
}

colnames(dff) <- c("i", "shannon")



plt <- dff %>% ggplot(aes(x = i, y = shannon)) +
    geom_line() +
    scale_x_continuous(breaks = round(seq(min(dff$i), max(dff$i), by = 500),1))

plt



pdf("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/secondary_str_shannon.pdf", width = 21, height = 7)
plt
dev.off()

