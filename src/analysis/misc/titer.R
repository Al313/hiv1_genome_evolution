

# load libraries
library(readxl)
library(stringr)
library(magrittr)
library(ggplot2)


## qpcr results
data <- read_excel("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/titration/qpcr.xls")
data <- data[!is.na(data$`Sample Name`),]

data$lines <- substring(data$SampleID, 8,9)
data$passage <- sapply(strsplit(data$SampleID, "p"), "[", 2)

data$Qunatity_per_2 <- data$`Quantity Mean` * 0.14 # to measure the copy number count of HIV-1 genomic RNA in 2


data$passage <- as.numeric(data$passage)
data$lines <- factor(data$lines, levels = c("13", "14", "15", "16"))

titration_fig <- data %>% ggplot(aes(x = passage, y = log10(Qunatity_per_2), color = lines)) +
                geom_line() +
                geom_point() +
                labs(title = "qPCR") +
                #ylim(c(1000000,40000000)) +
                scale_color_manual(values = c('red', 'blue', 'green', 'black')) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))


png(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/png/titration_qpcr.png")
titration_fig
dev.off()



## tcid


data <- as.data.frame(read_excel("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/titration/tcid50.xlsx"))


data$passage <- sapply(strsplit(data$Sample_ID, "p"), "[", 2)
data$passage[is.na(data$passage)] <- 0
data$lines <- substring(data$Sample_ID, 8,9)
data$lines[data$lines == ""] <- "ancestral"


data_sub <- data[data$Source_abv == "2022_rfc_i",]
data_sub <- data_sub[data_sub$passage != "NL4_3",]

data_sub$passage <- as.numeric(data_sub$passage)
data_sub$lines <- factor(data_sub$lines, levels = c("13", "14", "15", "16"))


tcid__fig1 <- data_sub %>% ggplot(aes(x = passage, y = log10(`TCID50/ml_MT2`), color = lines)) +
                geom_line() +
                geom_point() +
                #ylim(c(1000000,40000000)) +
                labs(title = "MT2_exp1") +
                scale_color_manual(values = c('red', 'blue', 'green', 'black')) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))




data_sub <- data[data$Source_abv == "2022_rfc_ii",]
data_sub <- data_sub[data_sub$passage != "NL4_3",]

data_sub$passage <- as.numeric(data_sub$passage)
data_sub$lines <- factor(data_sub$lines, levels = c("13", "14", "15", "16"))


tcid__fig2 <- data_sub %>% ggplot(aes(x = passage, y = log10(`TCID50/ml_MT2`), color = lines)) +
                geom_line() +
                geom_point() +
                #ylim(c(1000000,40000000)) +
                labs(title = "MT2_exp2") +
                scale_color_manual(values = c('red', 'blue', 'green', 'black')) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))






data_sub <- data[data$Source_abv == "2022_rfc_i",]
data_sub <- data_sub[data_sub$passage != "NL4_3",]

data_sub$passage <- as.numeric(data_sub$passage)
data_sub$lines <- factor(data_sub$lines, levels = c("13", "14", "15", "16"))


tcid__fig3 <- data_sub %>% ggplot(aes(x = passage, y = log10(`TCID50/ml_MT4`), color = lines)) +
                geom_line() +
                geom_point() +
                #ylim(c(1000000,40000000)) +
                labs(title = "MT4_exp1") +
                scale_color_manual(values = c('red', 'blue', 'green', 'black')) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))




data_sub <- data[data$Source_abv == "2022_rfc_ii",]
data_sub <- data_sub[data_sub$passage != "NL4_3",]

data_sub$passage <- as.numeric(data_sub$passage)
data_sub$lines <- factor(data_sub$lines, levels = c("13", "14", "15", "16"))


tcid__fig4 <- data_sub %>% ggplot(aes(x = passage, y = log10(`TCID50/ml_MT4`), color = lines)) +
                geom_line() +
                geom_point() +
                #ylim(c(1000000,40000000)) +
                labs(title = "MT4_exp2") +
                scale_color_manual(values = c('red', 'blue', 'green', 'black')) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))




pdf(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/tcid_100_500.pdf")


tcid__fig1
tcid__fig2
tcid__fig3
tcid__fig4

dev.off()




## moi


data_sub <- data[data$Source_abv %in% c("2022_rfc_i","2022_rfc_ii", "2018_christine","2024_christine"),]
data_sub <- data_sub[data_sub$passage != "NL4_3",]

data_sub$passage <- as.numeric(data_sub$passage)
data_sub$lines <- factor(data_sub$lines, levels = c("13", "14", "15", "16"))
data_sub$Source_abv <- factor(data_sub$Source_abv, levels = c("2024_christine","2018_christine","2022_rfc_i","2022_rfc_ii"))



moi__fig <- data_sub %>% ggplot(aes(x = passage, y = log10(MOI), color = lines)) +
                geom_line(aes(linetype=Source_abv)) +
                geom_point() +
                #ylim(c(1000000,40000000)) +
                labs(title = "moi") +
                scale_color_manual(values = c('red', 'blue', 'green', 'black')) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))


png(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/png/moi.png", height = 960, width = 960)
moi__fig
dev.off()
