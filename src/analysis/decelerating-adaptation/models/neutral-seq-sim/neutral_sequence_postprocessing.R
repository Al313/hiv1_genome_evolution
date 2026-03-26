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
total_run <- 10
long_df_all <- data.frame()
fixed_df_all <- data.frame(run_nr = integer(), exp_line = character(), passage = integer(), count = integer())

for (run_nr in 1:total_run){
  print(run_nr)
  fixed_df <- data.frame(run_nr = integer(), exp_line = character(), passage = integer(), count = integer())
  for (exp_line in exp_lines){
    print(exp_line)
    for (i in seq(1,50)){
      file_path_init <- paste0(wd, "results/tables/neutral-seq-sim/sequences/", bottleneck_freq, "/", exp_line, "/run", run_nr, "/init_population.npy")

      init_pop <- np$load(file_path_init)
      init_pop <- init_pop[1,]

      # print(i)
      file_path <- paste0(wd, "results/tables/neutral-seq-sim/sequences/", bottleneck_freq, "/", exp_line, "/run", run_nr, "/", i, ".npy")


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
      prop_mat_processed <- prop_mat[,apply(prop_mat, 2, function(col) any(col > 0.01)), drop = FALSE]

      prop_df_processed <- as.data.frame(prop_mat_processed)
      prop_df_processed[,"alt_allele"] <- rownames(prop_df_processed)

      if (ncol(prop_df_processed) == 1){
        next
      }

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
      long_df$run_nr <- run_nr
      long_df <- long_df[,c("run_nr","exp_line", "passage", "pos", "base_allele", "alt_allele", "allele_freq")]

      long_df_all <- rbind(long_df_all, long_df)
      count <- sum(long_df$allele_freq >= min_af)
      fixed_df <- rbind(fixed_df, data.frame(
		 			 run_nr = run_nr,
		  			 exp_line = exp_line,
		  			 passage = i*10,
		  			 count = count,
		  			 stringsAsFactors = FALSE))
    }
  }
  print(colnames(fixed_df))
  print(colnames(fixed_df_all))
  fixed_df_all <- rbind(fixed_df_all, fixed_df)
}

saveRDS(long_df_all, paste0(wd, "results/tables/neutral-seq-sim/", bottleneck_freq, "_sim_freq_all_runs.rds"))
saveRDS(fixed_df_all, paste0(wd, "results/tables/neutral-seq-sim/", bottleneck_freq, "_sim_count_all_runs.rds"))


