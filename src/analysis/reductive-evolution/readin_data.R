# loading libraries
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)

# determine the server path
if (file.exists("/home/amovas/")){
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
}

# define functions
`%notin%` <- Negate(`%in%`)

# load data
variants_ann_expiii <- read.table(file = paste0(wd, "results/tables/pipeline-outputs/iii/iii_annotated_variants.tsv.gz"), 
sep = "\t", stringsAsFactors = FALSE, header = TRUE)


# correct for the switched experimental lines
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == 15 & variants_ann_expiii$passage > 360] <- "A"
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == 16 & variants_ann_expiii$passage > 360] <- "B"

variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == "A"] <- "16"
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == "B"] <- "15"

variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == 13 & variants_ann_expiii$passage == 700] <- "A"
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == 16 & variants_ann_expiii$passage == 700] <- "B"

variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == "A"] <- "16"
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == "B"] <- "13"

# determine what data to be included
end_psg <- 600

# correct the feature name
variants_ann_expiii <- variants_ann_expiii %>% filter(passage <= end_psg) %>%
        mutate(feature = ifelse(feature == "envelope", "env", feature)) %>%
        mutate(feature = ifelse(feature == "3UTR", "U3", feature)) %>%
        mutate(feature = ifelse(feature == "5UTR", "U5", feature)) %>%
        mutate(feature = ifelse(feature == "5LTRLS", "GLS", feature))




# facors
host_factor <- c("MT2", "MT4")
exp_line_factor <- c("MT-2_1", "MT-2_2", "MT-4_1", "MT-4_2")
translational_impacts <- c("U", "S", "N")
gene_factor <- c("gag", "pol", "vif", "vpr", "tat", "rev", "vpu", "env", "nef")
feature_factor <- c("5R","U5", "GLS", "gag", "pol", "vif", "vpr", "tat1", "rev1", "vpu", "env", "tat2", "rev2", "nef", "U3", "3R")
feature_start <- c(454,551,635,790,2085,5041,5559,5830,5969,6061,6221,8369,8369,8787,9408,9529)
feature_end <- c(550,634,789,2292,5096,5619,5849,6044,6044,6306,8785,8414,8643,9407,9528,9625)

feature_df <- data.frame(feature = feature_factor,
                            start = feature_start,
                            end = feature_end,
                            length = feature_end - feature_start + 1)

feature_df_p <- rbind(feature_df, c("tat", sum(feature_df$length[feature_df$feature %in% c("tat1", "tat2")])))
feature_df_p$length <- as.numeric(feature_df_p$length)
feature_df_p <- rbind(feature_df_p, c("rev", sum(feature_df_p$length[feature_df_p$feature %in% c("rev1", "rev2")])))
feature_df_p$length <- as.numeric(feature_df_p$length)

feature_df_p <- feature_df_p[feature_df_p$feature %notin% c("tat1", "tat2", "rev1", "rev2"),]
feature_df_p$feature <- factor(feature_df_p$feature, levels = c("5R","U5", "GLS", "gag", "pol", "vif", "vpr", "tat", "rev", "vpu", "env", "nef", "U3", "3R"))



# format the dataframe
variants_ann_expiii <- variants_ann_expiii %>%
                            mutate(host_line = ifelse(exp_line %in% c("13", "14"), "MT2", "MT4"), # set host cell line info
                            exp_line = case_when(exp_line==13 ~ "MT-2_1", exp_line==14 ~ "MT-2_2", exp_line==15 ~ "MT-4_1", exp_line==16 ~ "MT-4_2"), # update exp line labels
                            mut_info = paste(ref_allele, genomic_pos, alt_allele, sep = "")) %>% 
                            mutate(effect_simplified = case_when(effect == "synonymous_variant" ~ "S", effect == "untranslated_region" ~ "U")) %>%
                            mutate(effect_simplified = replace_na(effect_simplified, "N")) %>%
                            mutate(passage = as.numeric(as.character(passage)),
                            exp_line = factor(exp_line, levels = exp_line_factor),
                            host_line = factor(host_line, levels = host_factor),
                            effect_simplified = factor(effect_simplified, levels = translational_impacts)) 


# remove annotation for U3 which overlaps with nef
variants_ann_expiii %<>% filter(!(genomic_pos >= 9076 & genomic_pos <= 9407 & effect_simplified == "U"))


# remove duplicate records for the same variants
variants_expiii <- variants_ann_expiii %>% filter(!duplicated(variant_id))
