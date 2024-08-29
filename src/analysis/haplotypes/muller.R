
# load libraries
library(ape)
library(stringr)
library(dplyr)
library(ggmuller)
library(mullerplot)

# read in MSA
msa <- read.dna("/Users/alimos313/Desktop/scrap/portal/all_msa.fasta", format = "fasta")

# calculate distance matrix
dist_matrix <- dist.dna(msa, model = "TN93")

# cluster
hc <- hclust(dist_matrix, method = "complete")

# save cluster tree
pdf(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/haplo_clustering.pdf", width = 100, height = 100)
plot(hc, main = "Haplotype Clustering Dendrogram", cex = 3)
dev.off()


# to view distance matrix for certain generations
mat <- as.matrix(dist_matrix)
mat[str_detect(row.names(mat), pattern = "460"),]



## first cutting threshold

height_threshold <- 0.0089

### OG

# define clades
clades <- cutree(hc, h = height_threshold)


# save cluster tree with clades
pdf(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/haplo_clustering_with_clades.pdf", width = 100, height = 100)
plot(hc, main = "Haplotype Clustering with Clades", cex = 3)
rect.hclust(hc, h = height_threshold, border = 2:4) 
dev.off()

# convert to data frame
haplotype_clades <- data.frame(
  Haplotype = names(clades),  # Row names from the original data
  Clade = clades
)

# modify class types
haplotype_clades$passage <- sapply(str_split(haplotype_clades$Haplotype, pattern = "_"), "[", 2)
haplotype_clades$abundance <- sapply(str_split(haplotype_clades$Haplotype, pattern = "_"), "[", 4)
haplotype_clades$abundance <- as.numeric(haplotype_clades$abundance)


# get abundance of clades per generation
haplotype_clades2 <- haplotype_clades %>% 
            group_by(Clade, passage) %>%
            summarize(cumulative_abundance = sum(abundance)) %>%
            ungroup()



# change colnames to be compatible with get_Muller_df function from ggmuller package
colnames(haplotype_clades2) <- c("Identity", "Generation", "Population")

# modify class types
haplotype_clades2$Generation <- as.integer(haplotype_clades2$Generation)
haplotype_clades2$Identity <- as.character(haplotype_clades2$Identity)



# define lineages and store in a compatible format
View(haplotype_clades2)
haplotype_clades2[haplotype_clades2$Generation == 230,]

edges <- data.frame(
  Parent = c("1", "2", "3", "3", "4", "6", "7"),
  Identity = c("2", "3", "4", "5", "6", "7", "8")
)



# get the muller data.frame and plot the Muller plot
muller_df <- get_Muller_df(edges, haplotype_clades2)
Muller_plot(muller_df, add_legend = TRUE)




pdf(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/haplo_muller.pdf", width = 14, height = 14)
Muller_plot(muller_df, add_legend = TRUE)
dev.off()




### with amendments



clades[clades==7] <- 6
clades[clades==8] <- 7

clades[clades==5] <- 4
clades[clades==6] <- 5
clades[clades==7] <- 6



edges <- data.frame(
  Parent = c("1", "2", "3", "4", "5"),
  Identity = c("2", "3", "4", "5", "6")
)


# tweaking for presentation purposes

haplotype_clades2$Population[haplotype_clades2$Generation == 470 & haplotype_clades2$Identity == 5] <- 0.82
haplotype_clades2$Population[haplotype_clades2$Generation == 470 & haplotype_clades2$Identity == 4] <- 0.18

# get the muller data.frame and plot the Muller plot
muller_df <- get_Muller_df(edges, haplotype_clades2)

plot <- Muller_plot(muller_df, add_legend = TRUE, conceal_edges = TRUE)




pdf(file = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/results/figs/pdf/haplo_muller_amended.pdf", width = 14, height = 14)
Muller_plot(muller_df, add_legend = TRUE)
dev.off()


## second cutting threshold


edges <- data.frame(
  Parent = c("1", "2", "3", "3", "4", "4", "5", "5", "5"),
  Identity = c("2", "3", "4", "5", "6", "7", "8", "9", "10")
)
