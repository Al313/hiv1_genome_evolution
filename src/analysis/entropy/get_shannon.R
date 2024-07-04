
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
    cl_line <- args[1]
}





if (dir.exists("/Users/alimos313")){
    main_wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/"
} else {
    main_wd <- "/home/amovas/data/"
}


# read in the variant data
variants_ann <- read.table(file = paste0(main_wd, "genome-evo-proj/results/tables/pipeline-outputs/all_annotated_variants.tsv.gz"), 
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




# define function to calculate Shannon entropy
shannon <- function(x, seq_err = 0.01){
    if(x >= seq_err) {-x*log2(x)} else (0)
}

# define function to determine polymorphism

polymorph <- function(x){
    ifelse(any(x >= 0.01), T, F)
}



# create a dataframe to store Shannon entropy results
diversity_df <- data.frame()
head(diversity_df)



line <- cl_line
variants_sub <- variants[variants$exp_line == line,]

for (psg in sort(unique(variants_sub$passage))){

    variants_sub2 <- variants_sub[variants_sub$passage == psg,]
    variants_sub2 <- variants_sub2[variants_sub2$genomic_pos >= 454 & variants_sub2$genomic_pos <= 470,]

    for (i in seq(454,470)){
        if (i %in% variants_sub2$genomic_pos){
            ss <- variants[variants$exp_line == line & variants$passage == psg & variants$genomic_pos == i,]
            alt_shannon <- sum(sapply(as.numeric(ss$allele_freq), FUN = shannon))
            ref_shannon <- shannon(1-sum(ss$allele_freq[ss$allele_freq >= 0.01]))
            total_shannon <- alt_shannon + ref_shannon

            polymorphic <- polymorph(ss$allele_freq)
        } else {
            total_shannon <- 0
            polymorphic <- 0
        }
        
        diversity_df <- rbind(diversity_df, c(line, psg, i, polymorphic, total_shannon))
    }

}



colnames(diversity_df) <- c("exp_line", "passage", "genomic_position", "polymorphic", "shannon_entropy")


write.table(diversity_df, file = paste0(main_wd, "genome-evo-proj/results/tables/misc/", cl_line,"_shannon_entropy.tsv"), quote = F, row.names = F, sep = "\t")


