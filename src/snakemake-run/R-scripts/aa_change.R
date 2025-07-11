
### load libraries

library(vcfR)
library(dplyr)
library(stringr)

### read-in the data

# determine the server path
if (file.exists("/home/amovas/")){
  print("Remote HPC Connection!")
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  print("Local PC Connection!")
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
}

# get cl arguments
args <- commandArgs(trailingOnly = TRUE)

# Access the first argument
exp <- args[1]

# read in cds feature list 

features_cds <- read.table(file = paste0(wd,"results/tables/pipeline-outputs/cds_feature_list.tsv"), 
                            sep = "\t", header = T, stringsAsFactors = F)


# read in annotated variants in vcf format to extract the INFO filed

ann_var_vcf <- read.vcfR(file = paste0(wd,"results/tables/pipeline-outputs/", exp, "/", exp, "_variants.ann.vcf.gz"))
ann_var_info <- vcfR::INFO2df(ann_var_vcf)
#ann_var_info <- ann_var_info[,-(6:7)]
#head(ann_var_info[ann_var_info$AF>= 0.9,], n = 1)


# read in the variant table

ann_var_df <- read.table(file = paste0(wd,"results/tables/pipeline-outputs/", exp, "/", exp, "_variants.ann.vcf.gz"), header = F,
                          sep = "\t", stringsAsFactors = F)
colnames(ann_var_df) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
ann_var_df <- ann_var_df[,1:5]

# merge the two dataset
ann_var <- cbind(ann_var_df, ann_var_info)



# add annotation of untranslated regions manually
ann_var$ANN[ann_var$POS >= 551 & ann_var$POS <= 634] <- paste0(ann_var$ANN[ann_var$POS >= 551 & ann_var$POS <= 634], ",", ann_var$ALT[ann_var$POS >= 551 & ann_var$POS <= 634], "|untranslated_region|MODERATE|5_UTR|5_UTR||AF324490.5|||||||||")
ann_var$ANN[ann_var$POS >= 9076 & ann_var$POS <= 9528] <- paste0(ann_var$ANN[ann_var$POS >= 9076 & ann_var$POS <= 9528], ",", ann_var$ALT[ann_var$POS >= 9076 & ann_var$POS <= 9528], "|untranslated_region|MODERATE|3_UTR|3_UTR||AF324499.3|||||||||")

# add annotation of additional untranslated regions manually
ann_var$ANN[ann_var$POS >= 454 & ann_var$POS <= 550] <- paste0(ann_var$ANN[ann_var$POS >= 454 & ann_var$POS <= 550], ",", ann_var$ALT[ann_var$POS >= 454 & ann_var$POS <= 550], "|untranslated_region|MODERATE|5_R|5_R||AF324489.5|||||||||")
ann_var$ANN[ann_var$POS >= 635 & ann_var$POS <= 789] <- paste0(ann_var$ANN[ann_var$POS >= 635 & ann_var$POS <= 789], ",", ann_var$ALT[ann_var$POS >= 635 & ann_var$POS <= 789], "|untranslated_region|MODERATE|5_LTR_LS|5_LTR_LS||AF324491.5|||||||||")
ann_var$ANN[ann_var$POS >= 9529 & ann_var$POS <= 9625] <- paste0(ann_var$ANN[ann_var$POS >= 9529 & ann_var$POS <= 9625], ",", ann_var$ALT[ann_var$POS >= 9529 & ann_var$POS <= 9625], "|untranslated_region|MODERATE|3_R|3_R||AF324500.3|||||||||")


### get the name of variant classes that you want to keep

# first unlist
ann <- unlist(str_split(ann_var$ANN, pattern = ","))

# get all the variant class category names
variant_calss <- data.frame(class = sapply(str_split(ann,pattern = "\\|"), "[[", 2), effect = sapply(str_split(ann,pattern = "\\|"), "[[", 3))

# get the counts in a df format
variant_calss_count <- as.data.frame(table(variant_calss$class,variant_calss$effect))
colnames(variant_calss_count) <- c("class", "effect", "count")

# filter out the ones that are insignificant
variant_calss_count <- variant_calss_count[variant_calss_count$effect != "MODIFIER" & variant_calss_count$count != 0,]

# store the category names that we would like to include in the annotation
included_variant_class <- as.vector(unique(variant_calss_count$class))



### extract the annotations of interest

ann <- str_split(ann_var$ANN, pattern = ",") #[c(1,85111,100001)]

# the dot is added so that for variants where there are multiple annotations of interest
# e.g.(the mutation occurs in a region where two genes overlap) when unlisting later on
# the number of annotation will be added after the dot 
#(and not just like a digit to the ID number which will not be discernible later whether this is part of the number or extra annotation)
names(ann) <- paste0(ann_var$ID,".") #[c(1,85111,100001)]

# make a filter boolean vector of the annotations that we want to include
ann_of_interest_bool <- lapply(ann, function(x) str_detect(x,str_c(included_variant_class, collapse = "|")))

# filter out the annotations that are not of interest
ann_of_interest <-  mapply(function(x,y) tryCatch(expr = x[y],error = {function(e) "ERROR"}),
                             x=ann, y=ann_of_interest_bool)

# remove the variants for which there is no annotation of interest
ann_of_interest <- ann_of_interest[!sapply(ann_of_interest, identical, character(0))]

### extract required information from the annotations

# unlsit the annotations to store them in a df (the added dot to the names is handy here)
ann_of_interest <- unlist(ann_of_interest, use.names = T)

# get the annotation fileds that are of interest
ann_df <- as.data.frame(do.call(rbind, lapply(sapply(ann_of_interest, str_split, "\\|"), "[",c(2,3,7,11,12,14))))

colnames(ann_df) <- c("effect", "impact", "feature_id", "aa_change", "loc_in_cds", "loc_in_protein")

# assign the variant and annotation id
ann_df$variant_id <- sapply(str_split(row.names(ann_df), pattern = "\\."), "[", 1)
ann_df$annotation_id <- row.names(ann_df)


row.names(ann_df) <- 1:nrow(ann_df)



### merging variants and annotations

# merge variant info with annotation info
ann_var_processed <- merge(ann_var[,-11], ann_df, by.x = "ID", by.y = "variant_id")


# merge variant and annotation info with gene info
ann_var_processed <- merge(ann_var_processed, features_cds[,c("feature_id","feature","ORF", "feature_range")], by="feature_id")

# reorganize the dataframe
colnames(ann_var_processed)[2:13] <- c("variant_id", "chrom", "genomic_pos", "ref_allele", "alt_allele", "allele_freq", "mut_type", "passage", "exp_line", "coverage", "loss_of_function", "nonsense_mediated_decay")
ann_var_processed <- ann_var_processed[,c(3:11,1,20,22,21,2,19,14:18,12:13)]
ann_var_processed <- ann_var_processed[order(ann_var_processed$genomic_pos,ann_var_processed$allele_freq),]


# change the chrom name to NL43_ann_wk0virusPassRef_plasmid
ann_var_processed[,1] <- "NL43_ann_wk0virusPassRef_plasmid"

# save the data frame
write.table(ann_var_processed, file = gzfile(paste0("/home/amovas/data/genome-evo-proj/results/tables/pipeline-outputs/", exp, "/", exp, "_annotated_variants.tsv.gz")),
    sep = "\t", quote = F, row.names = F)




### inspecting the data


# tail(ann_var_processed)

# min(ann_var_processed$genomic_pos)
# head(ann_var_processed[ann_var_processed$genomic_pos == 5580,])

# fixed <- ann_var_processed[ann_var_processed$allele_freq >= 0.99,]

# fixed[fixed$variant_id == 76680,]
# fixed[fixed$variant_id == 144103,]

# fixed[fixed$genomic_pos<= 8414 & fixed$genomic_pos >= 8369,]
