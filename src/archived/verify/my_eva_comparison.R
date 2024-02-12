


## in the following chunk it is shown that Eva's and my scripts produce exactly the same mutation calls but in different formats

eva_version_muts <- read.table(file = "/home/amovas/data/genome-evo-proj/results/tables/2-p/all_muts.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
eva_version_muts$unique_id <- paste(eva_version_muts$pos,eva_version_muts$new, eva_version_muts$line, eva_version_muts$passage, sep = "_")
eva_version_muts <- eva_version_muts[order(eva_version_muts$unique_id),]
#eva_version_muts <- eva_version_muts[eva_version_muts$line == 20 & eva_version_muts$passage %in% c(180,410),]


my_version_all <- read.table(file = "/home/amovas/data/genome-evo-proj/results/tables/3-p/all_variants.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
my_version_muts <- my_version_all[my_version_all$mut_type == "M",]
my_version_muts$unique_id <- paste(my_version_muts$start, my_version_muts$alt, my_version_muts$line, my_version_muts$passage, sep = "_")
my_version_muts <- my_version_muts[order(my_version_muts$unique_id),]

table(eva_version_muts$unique_id %in% my_version_muts$unique_id)


## let's compare the insertion results

eva_version_ins <- read.table(file = "/home/amovas/data/genome-evo-proj/results/tables/2-p/all_ins.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
eva_version_ins$unique_id <- paste(eva_version_ins$inserted, eva_version_ins$pos, eva_version_ins$line, sep = "_")
#eva_version_ins <- eva_version_ins[eva_version_ins$line == 20 & eva_version_ins$passage %in% c(180,410),]



my_version_ins <- my_version_all[my_version_all$mut_type == "I",]
my_version_ins$unique_id <- paste(my_version_ins$alt, my_version_ins$start, my_version_ins$line, sep = "_")

head(eva_version_ins)
head(my_version_ins)


table(eva_version_ins[eva_version_ins$fraction >= 0.99,"unique_id"])
nrow(eva_version_ins)
table(my_version_ins[my_version_ins$fraction >= 0.9,"unique_id"])
nrow(my_version_ins)

my_version_ins[my_version_ins$passage == 180 & my_version_ins$line == 20 & my_version_ins$fraction >= 0.8,]

## deletions


my_version_del <- my_version_all[my_version_all$mut_type == "D",]
my_version_del$unique_id <- paste(my_version_del$alt, my_version_del$start, my_version_del$line, sep = "_")

head(my_version_del)
my_version_del[my_version_del$fraction >= 0.8,]

my_version_del[my_version_del$passage == 180 & my_version_del$line == 20 & my_version_del$fraction >= 0.8,]
