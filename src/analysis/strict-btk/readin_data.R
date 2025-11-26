


# determine the server path

if (file.exists("/home/amovas/")){
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
}

# define functions

`%notin%` <- Negate(`%in%`)

# load exp-2e data
variants_ann_expiii <- read.table(file = paste0(wd, "results/tables/pipeline-outputs/iii/iii_annotated_variants.tsv.gz"), 
sep = "\t", stringsAsFactors = FALSE, header = TRUE)



# correct for the switched experimental lines
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == 13 & variants_ann_expiii$passage == 700] <- "A"
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == 15 & variants_ann_expiii$passage == 700] <- "B"

variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == "A"] <- "15"
variants_ann_expiii$exp_line[variants_ann_expiii$exp_line == "B"] <- "13"

# load exp-btk data
variants_ann_expiv <- read.table(file = paste0(wd, "results/tables/pipeline-outputs/iv/iv_annotated_variants.tsv.gz"), 
sep = "\t", stringsAsFactors = FALSE, header = TRUE)


# correct for the switched experimental lines
variants_ann_expiv$exp_line[variants_ann_expiv$exp_line == 18 & variants_ann_expiv$passage > 350] <- "A"
variants_ann_expiv$exp_line[variants_ann_expiv$exp_line == 19 & variants_ann_expiv$passage > 350] <- "B"

variants_ann_expiv$exp_line[variants_ann_expiv$exp_line == "A"] <- "19"
variants_ann_expiv$exp_line[variants_ann_expiv$exp_line == "B"] <- "18"

variants_ann_expiv <- variants_ann_expiv[!(variants_ann_expiv$exp_line == 20 & variants_ann_expiv$passage == 680),] # remove the 680th passage for the line 20, as it is low quality data


variants_ann_2exp <- rbind(variants_ann_expiii, variants_ann_expiv)

# determine what data to be included

end_psg <- 630

variants_ann_2exp <- variants_ann_2exp %>% filter(passage <= end_psg) %>%
        mutate(feature = ifelse(feature == "envelope", "env", feature)) # correct the feature name




# facors
exp_factor <- c("RB", "SB")
exp_line_factor <- c("MT-2_1", "MT-2_2", "MT-2_i", "MT-2_ii", "MT-2_iii", "MT-2_iv")
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
variants_ann_2exp <- variants_ann_2exp %>%
                            filter(exp_line %in% c(13:14,17:20) & allele_freq >= 0.01) %>%
                            mutate(exp_dataset = ifelse(exp_line %in% c("13", "14"), "RB", "SB"), # set host cell line info
                            exp_line = case_when(exp_line==13 ~ "MT-2_1", exp_line==14 ~ "MT-2_2", exp_line==17 ~ "MT-2_i", exp_line==18 ~ "MT-2_ii", exp_line==19 ~ "MT-2_iii", exp_line==20 ~ "MT-2_iv"), # update exp line labels
                            mut_info = paste(genomic_pos, ref_allele, alt_allele, sep = "_")) %>% 
                            mutate(effect_simplified = case_when(effect == "synonymous_variant" ~ "S", effect == "untranslated_region" ~ "U")) %>%
                            mutate(effect_simplified = replace_na(effect_simplified, "N")) %>%
                            mutate(exp_dataset = factor(exp_dataset, levels = exp_factor),
                              passage = as.factor(passage),
                              exp_line = factor(exp_line, levels = exp_line_factor),
                              effect_simplified = factor(effect_simplified, levels = translational_impacts)) 

# remove duplicate records for the same variants
variants_2exp <- variants_ann_2exp %>% filter(!duplicated(variant_id))






