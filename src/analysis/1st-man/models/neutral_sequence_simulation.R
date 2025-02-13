# get arguments from the commandline
args <- commandArgs(trailingOnly = TRUE)  
exp_line <- args[1]  
mut_cat <- args[2]
generation_time <- args[3]



# load libraries
library("stringr")


# determine the server path
if (file.exists("/home/amovas/")){
  print("Remote HPC Connection!")
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  print("Local PC Connection!")
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
  exp_line <- "MT-2_1"
  mut_cat <- "Majority"
  generation_time <- 180
}


# set parameters
category <- mut_cat
print(category)
if (category == "Sporadic"){
    min_threshold <- 0.01
    max_threshold <- 0.05
} else if (category == "Minority"){
    min_threshold <- 0.05
    max_threshold <- 0.5
} else if (category == "Majority"){
    min_threshold <- 0.5
    max_threshold <- 1.01
} else if (category == "Fixed"){
    min_threshold <- 0.99
    max_threshold <- 1.01
} else if (category == "All"){
    min_threshold <- 0.01
    max_threshold <- 1.01
}

## define functions

# define logarithmic mean function
logarithmic_mean <- function(a, b) {
  if (a == b) {
    return(a)  # The logarithmic mean of two equal numbers is just the number itself
  } else {
    return((a - b) / (log(a) - log(b)))
  }
}


# define transfer_size calculater function
calculate_transfer_size <- function(t, base_transfer_sizes) {
  if (t < 1 || t > 500) {
    stop("Time t must be between 1 and 500")
  }

  if (t <= 100) {
    interval = 10
  } else {
    interval = 100
  }

  floor_t = floor(t/interval) * interval
  ceil_t = ceiling(t/interval) * interval

  # Get transfer sizes for floor and ceiling time points
  size_floor = base_transfer_sizes[as.character(floor_t)]
  size_ceil = base_transfer_sizes[as.character(ceil_t)]

  result = logarithmic_mean(size_ceil,size_floor)

  return(round(as.numeric(result)))
}

# Function to mutate a base (returns integer)
modify_base <- function(x) {
  return(sample(setdiff(1:4, x), 1))  # Pick a different base
}


# use base encoding to save memory
base_encoding <- c(A = 1, C = 2, G = 3, T = 4)
base_decoding <- c("A", "C", "G", "T")

value_types <- as.factor(base_encoding)

# Function to calculate proportions while ensuring all categories are included
calculate_proportions <- function(col) {
  tab <- tabulate(factor(col, levels = value_types), nbins = length(value_types))
  tab / sum(tab)  # Convert to proportions
}


# Function to find values with frequency > 0.95 or return NA
find_dominant_values <- function(col, threshold = min_threshold) {
  dominant_values <- which(col > threshold)  # Get values with >99% proportion
  
  if (length(dominant_values) == 0) {
    return(NA)  # Return NA if no dominant values exist
  } else {
    return(dominant_values)  # Return the dominant value(s)
  }
}



# read-in MOI data for transfer size estimation

MOI_results <- readRDS(file = paste0(wd, "results/tables/moi_cleaned.rds"))

base_transfer_sizes <- MOI_results$btk_size[str_detect(MOI_results$exp_line, pattern = exp_line)][1:14]
names(base_transfer_sizes) <- unique(MOI_results$passage)
base_transfer_sizes <- append(base_transfer_sizes[1],base_transfer_sizes)
names(base_transfer_sizes)[1] <- "0"



# Set simulation parameters
host <- exp_line
print(host)
genome_length <- 9171   # HIV-1 genome length
print(genome_length)
initial_population <- 400  # Initial number of individuals
R0 <- 44  # Number of offspring per genome per generation
mutation_rate <- 2.16*(10^-5)  # Per base per replication
total_generations <- as.integer(generation_time)  # Total generations
print(total_generations)
bottleneck_intervals <- 2  # Every 2 generations, apply bottleneck





# Initialize population as an **integer matrix**
init_population <- matrix(rep(sample(1:4, genome_length, replace = T), initial_population), byrow = TRUE,
                          nrow = initial_population, ncol = genome_length)




# Store total mutations over time
mutation_counts <- numeric(round(total_generations/20))

set.seed(420)  # For reproducibility

# Simulation loop
for (gen in 1:total_generations) {
    print(gen)
    # set initial population
    if (gen == 1) {
        population <- init_population
    }

    bottleneck_size <- calculate_transfer_size(ceiling(gen/2), base_transfer_sizes)
    print(bottleneck_size)

    # Step 1: Mutation - Each base mutate with a probability of mutation_rate


    mutation_events <- which(runif(length(population)) < mutation_rate, arr.ind = TRUE)
    if (length(mutation_events) > 0) {
        population[mutation_events] <- sapply(population[mutation_events], modify_base)
    }

    # Step 2: Replication - Each genome produces R0 offspring

    population <- population[sample(nrow(population), nrow(population)*R0, replace = TRUE), ]  # Select parents
    


    # Step 3: Apply bottleneck every 2 generations
    if (gen %% bottleneck_intervals == 0) {
      population <- population[sample(nrow(population), bottleneck_size, replace = FALSE), ]
    }

    if (gen %% 20 == 0){
    psg <- round(gen/20)
    # Step 4: Determine variant frequency

    # Calculate proportions for each column
    proportions_list <- lapply(as.data.frame(population)[sample(1:nrow(population),round(nrow(population)/2)),], calculate_proportions)
    proportions <- do.call(cbind, proportions_list)

    # Apply function to each column of proportions_df
    dominant_values_per_column <- apply(proportions, 2, find_dominant_values)


    # Step 5: Track mutation accumulation
    mutation_counts[psg] <- length(which(as.character(dominant_values_per_column) != init_population[1,]))  # Total number of mutations across all genomes

    if (psg != 1){
            if (mutation_counts[psg] > mutation_counts[psg-1]){
                    print("######")
                    print(which(as.character(dominant_values_per_column) != init_population[1,]))
            }
    }
    }

}


write.table(mutation_counts, file = paste0(wd, "results/tables/misc/neutral-seq-sim/mut_counts/neutral_simulation_", host, "_", category, "_", genome_length, "_", total_generations, ".tsv"), quote = F, row.names = F, sep = "\t")
write.table(proportions, file = paste0(wd, "results/tables/misc/neutral-seq-sim/proportions/neutral_simulation_proportion", host, "_", genome_length, "_", total_generations, ".tsv"), quote = F, row.names = F, sep = "\t")


