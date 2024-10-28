
# define functions

`%notin%` <- Negate(`%in%`)

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


# for (i in 1:9173){
#     if (i %notin% df$i){
#         df <- rbind(df,c(i,median(df$shannon[df$i %in% max(1,i-5):min(9173,i+5)])))
#     }
# }



shannon_df <- data.frame()
for (i in 1:9173){
    shannon_df <- rbind(shannon_df, c(i, median(df$shannon[df$i %in% max(1,i-55):min(9173,i+55)], na.rm = T)))
    
}

colnames(shannon_df) <- c("i", "shannon")



plt <- shannon_df %>% ggplot(aes(x = i, y = shannon)) +
    geom_line() +
    scale_x_continuous(breaks = round(seq(min(shannon_df$i), max(shannon_df$i), by = 500),1))

plt



pdf("/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/secondary_str_shannon.pdf", width = 21, height = 7)
plt
dev.off()




### compare statistically

# read in the variant df


# ...



# compare with shape data

shape <- shape[,c("Nt", "med_55_nt")]
# align coordinates
shape$Nt <- shape$Nt + 454

plot(freq_df_all$freq, shape$med_55_nt)


# compare with bp data
# align coordinates
shannon_df$i <- shannon_df$i+454


variants_subset <- distinct(variants[variants$allele_freq >= 0.99,c("genomic_pos", "mut_info", "exp_line")])

freq_df <- as.data.frame(table(variants_subset$genomic_pos))
colnames(freq_df) <- c("i", "Freq")

freq_df$i <- as.numeric(as.character(freq_df$i))

freq_df_all <- data.frame()

for (i in 1:9627){
    if (i %notin% freq_df$i){
        freq_df_all <- rbind(freq_df_all,c(i,0))
    } else {
        freq_df_all <- rbind(freq_df_all, as.numeric(freq_df[freq_df$i == i,]))
    }
}

colnames(freq_df_all) <- c("i", "freq")
freq_df_all <- freq_df_all[-c(1:454),]

freq_df_all[freq_df_all$freq > 40,]
variants_subset[variants_subset$genomic_pos == 759,]

plot(freq_df_all$freq, shannon_df$shannon)





## categorize secondary structure data, count mutations in each region and perform fisher's exact test


# test shape

shape <- shape[,c("Nt", "med_55_nt")]
shape$Nt <- shape$Nt + 454


shape %<>% mutate(category_shape = case_when(med_55_nt <= as.numeric(quantile(shape$med_55_nt, probs = c(0.1,0.9))[1]) ~ "low_shape", 
                                        med_55_nt >= as.numeric(quantile(shape$med_55_nt, probs = c(0.1,0.9))[2]) ~ "high_shape")) %>%
            mutate(category_shape = ifelse(is.na(category_shape), "mid", category_shape))


variants_subset <-  variants %>% filter(host_line == "MT2", exp == "EXPIII", allele_freq >= 0.99) %>%
                    select(genomic_pos, mut_info) %>%
                    distinct()


cat_count <- table(shape$category)[1:2]
var_count <- table(shape$category[shape$Nt %in% variants_subset$genomic_pos])[1:2]


fisher.test(cbind(cat_count,var_count))


# test shannon 

shannon_df$i <- shannon_df$i+454
shannon_df %<>% mutate(category_shannon = case_when(shannon <= as.numeric(quantile(shannon_df$shannon, probs = c(0.1,0.9))[1]) ~ "low_shannon", 
                                        shannon >= as.numeric(quantile(shannon_df$shannon, probs = c(0.1,0.9))[2]) ~ "high_shannon")) %>%
            mutate(category_shannon = ifelse(is.na(category_shannon), "mid", category_shannon))

variants_subset <-  variants %>% filter(host_line == "MT2", exp == "EXPIII", allele_freq >= 0.05) %>%
                    select(genomic_pos, mut_info) %>%
                    distinct()





cat_count <- table(shannon_df$category)[1:2]
var_count <- table(shannon_df$category[shannon_df$i %in% variants_subset$genomic_pos])[1:2]

fisher.test(cbind(cat_count,var_count))



# merge shape and base-pairing shannon data


merged_df <- merge(shape, shannon_df, by.x = "Nt", by.y = "i")


variants_subset <-  variants %>% filter(host_line == "MT4", exp == "EXPIII", allele_freq >= 0.05) %>%
                    select(genomic_pos, mut_info) %>%
                    distinct()



cat_count <- table(merged_df$category_shape, merged_df$category_shannon)
var_count <- table(merged_df$category_shape[merged_df$Nt %in% variants_subset$genomic_pos],merged_df$category_shannon[merged_df$Nt %in% variants_subset$genomic_pos])


## compare between double-high and double-low regions
fisher.test(matrix(c(var_count[1,1], var_count[2,2], cat_count[1,1], cat_count[2,2]), ncol = 2, byrow = F))

## compare variant count between double-low region and the rest of genome
mat <- matrix(c(sum(var_count)-var_count[2,2], var_count[2,2], sum(cat_count)-cat_count[2,2], cat_count[2,2]), ncol = 2, byrow = F)
fisher.test(mat)

chisq.test(mat)$expected



### are there bias in transition:transversion ratio?

# compute for whole genome

variants_tstv <-    variants %>% filter(mut_type == "P") %>%
                    mutate(tstv = case_when(ref_allele == "A" & alt_allele == "G" ~ "ts",
                                            ref_allele == "G" & alt_allele == "A" ~ "ts",
                                            ref_allele == "C" & alt_allele == "T" ~ "ts",
                                            ref_allele == "T" & alt_allele == "C" ~ "ts")) %>%
                    mutate(tstv = ifelse(is.na(tstv), "tv", tstv))

variants_subset <-  variants_tstv %>% filter(host_line == "MT2", exp == "EXPIII", allele_freq >= 0.05) %>%
                    select(genomic_pos, mut_info, tstv) %>%
                    distinct()


head(variants_subset)
table(variants_subset$tstv[variants_subset$genomic_pos %notin% c(455:789, 9076:9625)])

var_count <- table(merged_df$category_shape[merged_df$Nt %in% variants_subset$genomic_pos[variants_subset$tstv == "ts" & variants_subset$genomic_pos %in% c(455:789, 9076:9625)]],merged_df$category_shannon[merged_df$Nt %in% variants_subset$genomic_pos[variants_subset$tstv == "ts" & variants_subset$genomic_pos %in% c(455:789, 9076:9625)]])
var_count <- table(merged_df$category_shape[merged_df$Nt %in% variants_subset$genomic_pos[variants_subset$tstv == "tv" & variants_subset$genomic_pos %in% c(455:789, 9076:9625)]],merged_df$category_shannon[merged_df$Nt %in% variants_subset$genomic_pos[variants_subset$tstv == "tv" & variants_subset$genomic_pos %in% c(455:789, 9076:9625)]])

merged_df$Nt[merged_df$category_shape == "low_shape" & merged_df$category_shannon == "low_shannon"]

genomic_pos %in% c(455:789, 9076:9625)




# in non-coding regions tstv between double-low and WG excluding double-low
fisher.test(matrix(c(174,92,13,9), byrow = F, ncol = 2))

# in coding regions tstv between double-low and WG excluding double-low
fisher.test(matrix(c(1609,387,25,7), byrow = F, ncol = 2))


