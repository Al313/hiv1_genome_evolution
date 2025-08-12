
# load libraries


library("magrittr")
library("stringr")
library("dplyr")
library("tidyr")


# take command line variables
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0){
    cl_line <- 13
} else {
    cl_line <- as.integer(args[1])
}

if (cl_line %in% 13:16){
	exp <- "iii"
} else {
	exp <- "iv"
}

print(cl_line)



if (dir.exists("/Users/alimos313")){
    main_wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/"
} else {
    main_wd <- "/home/amovas/data/"
}


# read in the variant data
variants_ann <- read.table(file = paste0(main_wd, "genome-evo-proj/results/tables/pipeline-outputs/", exp, "/", exp, "_annotated_variants.tsv.gz"), 
sep = "\t", stringsAsFactors = FALSE, header = TRUE)


# account for the switched experimental lines
variants_ann$exp_line[variants_ann$exp_line == 15 & variants_ann$passage > 360] <- "A"
variants_ann$exp_line[variants_ann$exp_line == 16 & variants_ann$passage > 360] <- "B"

variants_ann$exp_line[variants_ann$exp_line == "A"] <- "16"
variants_ann$exp_line[variants_ann$exp_line == "B"] <- "15"



variants_ann$exp_line[variants_ann$exp_line == 18 & variants_ann$passage > 350] <- "A"
variants_ann$exp_line[variants_ann$exp_line == 19 & variants_ann$passage > 350] <- "B"

variants_ann$exp_line[variants_ann$exp_line == "A"] <- "19"
variants_ann$exp_line[variants_ann$exp_line == "B"] <- "18"


## prepare data

variants_ann %<>% 
                filter(coverage >= 500) %>%
                filter(mut_type == "P") %>%
                mutate(exp = case_when(exp_line %in% 13:16 ~ "EXPIII", exp_line %in% 17:20 ~ "EXPIV")) %>% # add experiment number info
                mutate(host_line = ifelse(exp_line %in% c("13", "14","17", "18", "19", "20"), "MT2", "MT4")) %>% # add a field containing pos, ref, and alt info of muts
                mutate(mut_info = paste(genomic_pos, ref_allele, alt_allele, sep = "_")) %>% # add host cell line info
                mutate(mut_cat = case_when((allele_freq < 0.01) ~ "Unreliable",
                                        (allele_freq >= 0.01 & allele_freq < 0.05) ~ "Sporadic", 
                                        (allele_freq >= 0.05 & allele_freq < 0.5) ~ "Minority",
                                        (allele_freq >= 0.5 & allele_freq < 0.99) ~ "Majority", 
                                        (allele_freq >= 0.99) ~ "Fixated"))

variants <- variants_ann %>% filter(!duplicated(variant_id))




# vectorized Shannon entropy
shannon <- function(x, seq_err = 0.01) {
    ifelse(x >= seq_err, -x * log2(x), 0)
}

# vectorized heterozygosity helper
square <- function(x, seq_err = 0.01) {
    ifelse(x >= seq_err, x^2, 0)
}

# define function to determine polymorphism
polymorph <- function(x){
    ifelse(any(x >= 0.01), T, F)
}


# subset data
variants_sub <- variants[variants$exp_line == cl_line & variants$genomic_pos >= 454 & variants$genomic_pos <= 9625,]


# Compute metrics only for positions that exist
diversity_df <- variants_sub %>%
    filter(allele_freq >= 0.01) %>%
    group_by(exp_line, passage, genomic_pos) %>%
    summarise(
        alt_shannon = sum(shannon(as.numeric(allele_freq))),
        ref_shannon = shannon(1 - sum(allele_freq)),
        alt_hetero  = sum(square(as.numeric(allele_freq))),
        ref_hetero  = square(1 - sum(allele_freq)),
        polymorphic = polymorph(allele_freq),
        .groups = "drop"
    ) %>%
    mutate(
        shannon_entropy = alt_shannon + ref_shannon,
        heterozygosity  = 1 - ref_hetero - alt_hetero
    ) %>%
    select(exp_line, passage, genomic_position = genomic_pos, polymorphic, heterozygosity, shannon_entropy)



# Add missing positions back in
all_pos <- expand.grid(
    exp_line = as.character(unique(variants_sub$exp_line)),
    passage = sort(unique(variants_sub$passage)),
    genomic_position = 454:9625
)


diversity_df <- all_pos %>%
    left_join(diversity_df, by = c("exp_line", "passage", "genomic_position")) %>%
    mutate(
        polymorphic = replace_na(polymorphic, FALSE),
        heterozygosity = replace_na(heterozygosity, 0),
        shannon_entropy = replace_na(shannon_entropy, 0)
    ) %>%
    arrange(passage, genomic_position)





write.table(diversity_df, file = paste0(main_wd, "genome-evo-proj/results/tables/diversity/", exp, "/", cl_line,"_diversity_measures.tsv"), quote = F, row.names = F, sep = "\t")


