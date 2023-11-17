
# loading libraries
library("plotly")
library("ggplot2")
library("DT")
library("dplyr")
library("tidyr")
library("Biostrings")
library("stringr")
library("magrittr")
library("ggvenn")
library("VennDiagram")
library("ggVennDiagram")




# read in the data
variants_ann <- read.table(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/tables/2-p/all_annotated_variants.tsv.gz", 
sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# set correct data types
variants_ann$exp_line <- as.integer(variants_ann$exp_line)


# set host cell line info
variants_ann$host_line <- ifelse(variants_ann$exp_line %in% c("13", "14","17", "18", "19", "20"), "MT2", "MT4")


# account for the switched experimental lines
variants_ann$exp_line[variants_ann$exp_line == 15 & variants_ann$passage > 360] <- "A"
variants_ann$exp_line[variants_ann$exp_line == 16 & variants_ann$passage > 360] <- "B"

variants_ann$exp_line[variants_ann$exp_line == "A"] <- "16"
variants_ann$exp_line[variants_ann$exp_line == "B"] <- "15"



variants_ann$exp_line[variants_ann$exp_line == 18 & variants_ann$passage > 350] <- "A"
variants_ann$exp_line[variants_ann$exp_line == 19 & variants_ann$passage > 350] <- "B"

variants_ann$exp_line[variants_ann$exp_line == "A"] <- "19"
variants_ann$exp_line[variants_ann$exp_line == "B"] <- "18"


# add a field containing pos, ref, and alt info of muts
variants_ann$mut_info <- paste(variants_ann$genomic_pos, variants_ann$ref_allele, variants_ann$alt_allele, sep = "_")

## EXPIII data
threshold <- 0.5

variants_ann_expiii <- variants_ann[variants_ann$exp_line %in% 13:16,]
variants_ann_expiii <- variants_ann_expiii %>%
                        mutate(exp_line = case_when(exp_line==13 ~ "MT2_1", exp_line==14 ~ "MT2_2", exp_line==15 ~ "MT4_1", exp_line==16 ~ "MT4_2")) %>%
                        filter(allele_freq >= threshold)

variants_expiii <- variants_ann_expiii %>% filter(!duplicated(variant_id))



## EXPIV data

variants_ann_expiv <- variants_ann[variants_ann$exp_line %in% 17:20,]
variants_ann_expiv <- variants_ann_expiv %>%
                        mutate(exp_line = case_when(exp_line==17 ~ "MT2_1", exp_line==18 ~ "MT2_2", exp_line==19 ~ "MT2_3", exp_line==20 ~ "MT2_4")) %>%
                        filter(allele_freq >= threshold)

variants_expiv <- variants_ann_expiv %>% filter(!duplicated(variant_id))


end_passage = 560


# uniqueness of majority mutations present at the end passage
end_maj <- variants_expiii %>%
                        filter(passage == end_passage) %>%
                        select(c(exp_line,mut_info)) %>%
                        distinct()



#table(as.numeric(table(par$mut_info)))


# create a named vector from the df
muts <- end_maj$mut_info
names(muts) <- end_maj$exp_line

# convert named vector to list
muts_list <- lapply(split(muts, names(muts)), unname)


head(muts_list)

all_muts <- as.character(unlist(muts_list))


k <- 0
shared_df <- data.frame()

while (k < 101){
    k <- k + 1

    all_muts_randomized <- sample(all_muts)


    .length = as.numeric(lapply(muts_list, length))

    muts_list_randomized <- split(all_muts_randomized, rep.int(seq_along(.length), .length))
    names(muts_list_randomized) <- names(muts_list)

    muts_list_randomized <- lapply(muts_list_randomized, unique)






    shared_percent <- c()

    for (i in 1:4){



            dd <- muts_list_randomized[[i]]
            for (j in 1:4){
                    shared_percent <- append(shared_percent,100*(sum(dd %in% muts_list_randomized[[j]])/(length(muts_list_randomized[[i]])+length(muts_list_randomized[[j]])-sum(dd %in% muts_list_randomized[[j]]))))
            }
            
    }


    shared_mat <- matrix(shared_percent, ncol = 4)
    shared_df <- as.data.frame(shared_mat) %>% 
                            `colnames<-`(names(muts_list_randomized)) %>%
                            mutate(line1 = names(muts_list_randomized)) %>%
                            gather(key = "line2", value = "val", 1:4) %>%
                            mutate(val = round(val, 2))
    colnames(shared_df)[3] <- k
    if (k == 1){
        shared_df_all <- shared_df
    } else {
        shared_df_all <- cbind(shared_df_all, shared_df[,3])
    }

}

shared_df_all <- rowMeans(shared_df_all[,3:ncol(shared_df_all)])

shared_df_all

shared_percent <- c()

for (i in 1:4){



        dd <- muts_list[[i]]
        for (j in 1:4){
                shared_percent <- append(shared_percent,100*(sum(dd %in% muts_list[[j]])/(length(muts_list[[i]])+length(muts_list[[j]])-sum(dd %in% muts_list[[j]]))))
        }
        
}

shared_percent


shared_mat <- matrix(shared_percent, ncol = 4)
shared_df <- as.data.frame(shared_mat) %>% 
                        `colnames<-`(names(muts_list)) %>%
                        mutate(line1 = names(muts_list)) %>%
                        gather(key = "line2", value = "val", 1:4) %>%
                        mutate(val = round(val, 2))
