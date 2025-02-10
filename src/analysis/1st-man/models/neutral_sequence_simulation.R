

# set parameters
category <- "Fixed"

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

  size_t = base_transfer_sizes[as.character(t)]

  result = logarithmic_mean(size_ceil,size_floor)

  return(result)
}

# define base_transfer_sizes from MOI data

base_transfer_sizes <- NA

# determine the server path

if (file.exists("/home/amovas/")){
  print("Remote HPC Connection!")
  wd <- "/home/amovas/data/genome-evo-proj/"
} else {
  print("Local PC Connection!")
  wd <- "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
}

set.seed(42)  # For reproducibility

# Parameters
genome_length <- 9171   # HIV-1 genome length
initial_population <- 400  # Initial number of individuals
R0 <- 44  # Number of offspring per genome per generation
mutation_rate <- 2.16*(10^-5)  # Per base per replication
total_generations <- 1000  # Total generations
bottleneck_intervals <- 2  # Every 2 generations, apply bottleneck
# bottleneck_size <- 400  # Approximate number of viruses transferred

# Initialize population (matrix of genomes, each row is an individual)
init_population <- matrix(rep(sample(c("A","C","G","T"), genome_length, replace = T), initial_population), nrow = initial_population, ncol = genome_length, byrow = T)  # Start with no mutations

# Store total mutations over time
mutation_counts <- numeric(total_generations)

# Function to modify a base
modify_base <- function(x) {
    return(sample(setdiff(c("A", "C", "G", "T"), x),1))  # Modify value based on existing value
}

# Define possible values
value_types <- c("A", "C", "G", "T")

# Function to calculate proportions while ensuring all categories are included
calculate_proportions <- function(col) {
  tab <- table(factor(col, levels = value_types))  # Ensure all values exist
  prop.table(tab)  # Convert to proportions
}


# Function to find values with frequency > 0.95 or return NA
find_dominant_values <- function(col) {
  dominant_values <- names(col[col > min_threshold])  # Get values with >99% proportion
  
  if (length(dominant_values) == 0) {
    return(NA)  # Return NA if no dominant values exist
  } else {
    return(dominant_values)  # Return the dominant value(s)
  }
}

# Simulation loop
for (gen in 1:total_generations) {
    print(gen)
    # set initial population
    if (gen == 1) {
        population <- init_population
    }

    bottleneck_size <- calculate_transfer_size(gen, base_transfer_sizes)

    # Step 2: Mutation - Each base mutates with probability mutation_rate
    mutations <- matrix(runif(nrow(population) * genome_length) < mutation_rate, 
                        nrow = nrow(population), ncol = genome_length)
    population[mutations] <- sapply(population[mutations], modify_base)
    
    # Step 1: Replication - Each genome produces R0 offspring
    new_population_size <- nrow(population) * R0
    new_population <- population[sample(nrow(population), new_population_size, replace = TRUE), ]  # Select parents
    
    
    # Step 3: Apply bottleneck every 2 generations
    if (gen %% bottleneck_intervals == 0) {
        new_population <- new_population[sample(nrow(new_population), bottleneck_size, replace = FALSE), ]
    } 
    
    # Step 4: Update population
    population <- new_population
    
    # Step 5: Track mutation accumulation


    # Calculate proportions for each column
    proportions <- apply(population, 2, calculate_proportions)

    # Convert to a data frame for better readability
    proportions_df <- as.data.frame(proportions)


    # Apply function to each column of proportions_df
    dominant_values_per_column <- apply(proportions_df, 2, find_dominant_values)
    # proportions_df[,is.na(as.character(dominant_values_per_column) == init_population[1,])]

    mutation_counts[gen] <- length(which(as.character(dominant_values_per_column) != init_population[1,]))  # Total number of mutations across all genomes
}


saveRDS(mutation_counts, file = paste0(wd, "results/tables/misc/neutral-seq-sim/neutral_simulation_", category, "_", genome_length, ".rds"))



