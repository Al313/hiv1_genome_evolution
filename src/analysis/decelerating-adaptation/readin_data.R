# loading libraries
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)



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

# determine what data to be included

end_psg <- 500

variants_ann_expiii <- variants_ann_expiii %>% filter(passage <= end_psg) %>%
        mutate(feature = ifelse(feature == "envelope", "env", feature)) # correct the feature name




# facors
host_factor <- c("MT2", "MT4")
exp_line_factor <- c("MT-2_1", "MT-2_2", "MT-4_1", "MT-4_2")
translational_impacts <- c("U", "S", "N")
gene_factor <- c("gag", "pol", "vif", "vpr", "tat", "rev", "vpu", "env", "nef")
feature_factor <- c("5R","5UTR", "5LTRLS", "gag", "pol", "vif", "vpr", "tat1", "rev1", "vpu", "env", "tat2", "rev2", "nef", "3UTR", "3R")
feature_start <- c(454,551,635,790,2085,5041,5559,5830,5969,6061,6221,8369,8369,8787,9408,9529)
feature_end <- c(550,634,789,2292,5096,5619,5849,6044,6044,6306,8785,8414,8643,9407,9528,9625)

feature_df <- data.frame(feature = feature_factor,
                            start = feature_start,
                            end = feature_end,
                            length = feature_end - feature_start + 1)

feature_df <- rbind(feature_df, c("tat", sum(feature_df$length[feature_df$feature %in% c("tat1", "tat2")])))
feature_df$length <- as.numeric(feature_df$length)
feature_df <- rbind(feature_df, c("rev", sum(feature_df$length[feature_df$feature %in% c("rev1", "rev2")])))
feature_df$length <- as.numeric(feature_df$length)

feature_df <- feature_df[feature_df$feature %notin% c("tat1", "tat2", "rev1", "rev2"),]
feature_df$feature <- factor(feature_df$feature, levels = c("5R","5UTR", "5LTRLS", "gag", "pol", "vif", "vpr", "tat", "rev", "vpu", "env", "nef", "3UTR", "3R"))



# format the dataframe
variants_ann_expiii <- variants_ann_expiii %>%
                            mutate(host_line = ifelse(exp_line %in% c("13", "14"), "MT2", "MT4"), # set host cell line info
                            exp_line = case_when(exp_line==13 ~ "MT-2_1", exp_line==14 ~ "MT-2_2", exp_line==15 ~ "MT-4_1", exp_line==16 ~ "MT-4_2"), # update exp line labels
                            mut_info = paste(ref_allele, genomic_pos, alt_allele, sep = "")) %>% 
                            mutate(effect_simplified = case_when(effect == "synonymous_variant" ~ "S", effect == "untranslated_region" ~ "U")) %>%
                            mutate(effect_simplified = replace_na(effect_simplified, "N")) %>%
                            mutate(passage = as.factor(passage),
                            exp_line = factor(exp_line, levels = exp_line_factor),
                            host_line = factor(host_line, levels = host_factor),
                            effect_simplified = factor(effect_simplified, levels = translational_impacts)) 

# remove duplicate records for the same variants
variants_expiii <- variants_ann_expiii %>% filter(!duplicated(variant_id))


line_col_palette <- c("#ff00ff", "#ff2400", "#6600cc", "#0000ff")
names(line_col_palette) <- c("MT-2_1", "MT-2_2", "MT-4_1", "MT-4_2")
