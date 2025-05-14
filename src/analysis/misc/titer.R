

# load libraries
library(readxl)
library(stringr)
library(magrittr)
library(ggplot2)



# qpcr expiii

line_col_palette <- c("#ff00ff", "#ff2400", "#6600cc", "#0000ff")




data <- read_excel("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/titration/SummaryMasterfile_qPCR_Titers.xlsx", sheet = 1, skip = 1)


data %<>% filter(exp_line %in% 13:16)

# correct for the switched experimental lines
data$exp_line[data$exp_line == 15 & data$passage > 360] <- "A"
data$exp_line[data$exp_line == 16 & data$passage > 360] <- "B"

data$exp_line[data$exp_line == "A"] <- "16"
data$exp_line[data$exp_line == "B"] <- "15"




# data[data$transferred_quantity_EXP_III <= 1000000,]



data$passage <- as.numeric(data$passage)
data$exp_line <- factor(data$exp_line, levels = c("13", "14", "15", "16", "17", "18", "19", "20"))


titration_fig <- data %>% ggplot(aes(x = passage, y = log10(transferred_quantity_EXP_III), color = exp_line)) +
                geom_line() +
                geom_point() +
                labs(title = "qPCR") +
                #ylim(c(1000000,40000000)) +
                scale_color_manual(values = line_col_palette) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))


png(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/png/titration_qpcr_expiii.png", height = 960, width = 960)
titration_fig
dev.off()

# qpcr expiii & expiv



line_col_palette <- c("#ff00ff", "#ff2400", "#6600cc", "#0000ff", "#66ff66", "#00cc99", "#009900", "#558000")




data <- read_excel("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/data/titration/SummaryMasterfile_qPCR_Titers.xlsx", sheet = 1, skip = 1)


data %<>% filter(exp_line %in% 13:20)

# correct for the switched experimental lines
data$exp_line[data$exp_line == 15 & data$passage > 360] <- "A"
data$exp_line[data$exp_line == 16 & data$passage > 360] <- "B"

data$exp_line[data$exp_line == "A"] <- "16"
data$exp_line[data$exp_line == "B"] <- "15"


data$exp_line[data$exp_line == 18 & data$passage > 350] <- "A"
data$exp_line[data$exp_line == 19 & data$passage > 350] <- "B"

data$exp_line[data$exp_line == "A"] <- "19"
data$exp_line[data$exp_line == "B"] <- "18"




data$passage <- as.numeric(data$passage)
data$exp_line <- factor(data$exp_line, levels = c("13", "14", "15", "16", "17", "18", "19", "20"))


titration_fig <- data %>% ggplot(aes(x = passage, y = log10(mean_quantity_per_ml_cell_culture_supernatant), color = exp_line)) +
                geom_line() +
                geom_point() +
                labs(title = "qPCR") +
                #ylim(c(1000000,40000000)) +
                scale_color_manual(values = line_col_palette) +
                theme_bw() +
                theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
                axis.title.x = element_text(size = 25),
                axis.title.y = element_text(size = 25),
                axis.text = element_text(size = 15))


png(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/png/titration_qpcr_expiii_iv.png", height = 960, width = 960)
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
