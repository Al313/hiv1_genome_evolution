#!/usr/bin/env Rscript

library(reticulate)
library(tidyr)

# Use Python's numpy to load the file
np <- import("numpy")

# determine the server path
if (file.exists("/home/amovas/")){
  print("Remote HPC Connection!")
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  print("Local PC Connection!")
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
}


exp_lines <- c("MT-2_1", "MT-2_2", "MT-4_1", "MT-4_2")
bottleneck_freq <- 3
min_af <- 0.99
mut_count_all <- list()
long_df_all <- data.frame()

for (exp_line in exp_lines){
  print(exp_line)
  mut_count <- c()
  for (i in seq(1,50)){
    file_path_init <- paste0(wd, "results/tables/neutral-seq-sim/sequences/", bottleneck_freq, "/", exp_line, "/init_population.npy")

    init_pop <- np$load(file_path_init)
    init_pop <- init_pop[1,]

    # print(i)
    file_path <- paste0(wd, "results/tables/neutral-seq-sim/sequences/", bottleneck_freq, "/", exp_line, "/", i, ".npy")


    pop <- np$load(file_path)


    # Define the possible values
    values <- 1:4

    # Initialize a matrix to store proportions
    prop_mat <- matrix(0, nrow = length(values), ncol = ncol(pop),
                      dimnames = list(values, 1:ncol(pop)))

    # Calculate proportions for each column
    for (j in seq_len(ncol(pop))) {
      # Count occurrences of each value in the column
      counts <- table(factor(pop[, j], levels = values))
      # Convert counts to proportions
      prop_mat[, j] <- counts / sum(counts)
    }


    # set ref allele freq to 0
    for (j in seq_along(init_pop)) {
      prop_mat[as.character(init_pop[j]), j] <- 0
    }


    # View the proportion matrix
    prop_mat_processed <- prop_mat[,apply(prop_mat, 2, function(col) any(col > 0.01))]

    prop_df_processed <- as.data.frame(prop_mat_processed)
    prop_df_processed[,"alt_allele"] <- rownames(prop_df_processed)


    long_df <- pivot_longer(
      data = prop_df_processed,
      cols = -starts_with("alt_allele"),  # Adjust this to match your measurement columns
      names_to = "pos",            # New column to indicate the original column names
      values_to = "allele_freq"           # New column to hold the values
    )

    long_df$pos <- as.integer(long_df$pos)
    long_df <- long_df[long_df$allele_freq != 0,]
    long_df$base_allele <- init_pop[long_df$pos]

    long_df <- long_df[order(long_df$pos),]
    long_df$passage <- i*10
    long_df$exp_line <- exp_line
    long_df <- long_df[,c("exp_line", "passage", "pos", "base_allele", "alt_allele", "allele_freq")]

    long_df_all <- rbind(long_df_all, long_df)
    count <- sum(long_df$allele_freq >= min_af)
    mut_count <- append(mut_count, count)
    mut_count_all[[exp_line]] <- mut_count
  }
}


saveRDS(long_df_all, paste0(wd, "results/tables/neutral-seq-sim/", bottleneck_freq, "_sim_freq_all1.rds"))
saveRDS(mut_count_all, paste0(wd, "results/tables/neutral-seq-sim/", bottleneck_freq, "_sim_count_all1.rds"))


