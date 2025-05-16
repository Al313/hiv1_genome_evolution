#!/usr/bin/env python3




import numpy as np
import pandas as pd
import sys
import os
import math

# Get command-line arguments
args = sys.argv[1:] if len(sys.argv) > 1 else ["MT-2_1", "Majority", 180]
exp_line, mut_cat, generation_time, sample_nr, bottleneck_freq, seq_sampling_freq = args

# Determine working directory
if os.path.exists("/home/amovas/"):
    print("Remote HPC Connection!")
    wd = "/home/amovas/data/genome-evo-proj/"
else:
    print("Local PC Connection!")
    wd = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
    exp_line, mut_cat, generation_time = "MT-4_2", "Majority", 180

# Set mutation category thresholds
category = mut_cat
print(category, flush = True)
category_thresholds = {
    "Sporadic": (0.01, 0.05),
    "Minority": (0.05, 0.5),
    "Majority": (0.5, 1.01),
    "Fixed": (0.99, 1.01),
    "All": (0.01, 1.01),
}

min_threshold, max_threshold = category_thresholds.get(category, (0.01, 1.01))

# Define logarithmic mean function
def logarithmic_mean(a, b):
    return a if a == b else (a - b) / (np.log(a) - np.log(b))

# Define transfer size calculation function
def calculate_transfer_size(t, base_transfer_sizes):
    if not (1 <= t <= 500):
        raise ValueError("Time t must be between 1 and 500")
    interval = 10 if t < 100 else 100
    floor_t, ceil_t = math.floor(t / interval) * interval, math.ceil(t / interval) * interval
    size_floor, size_ceil = base_transfer_sizes.get(floor_t), base_transfer_sizes.get(ceil_t)
    return round(logarithmic_mean(size_ceil, size_floor))

# Function to mutate a base
def modify_base(x):
    return np.random.choice([b for b in range(1, 5) if b != x])

# Base encoding to save memory
base_encoding = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
base_decoding = {v: k for k, v in base_encoding.items()}
value_types = np.array(list(base_encoding.values()))

# Function to calculate proportions ensuring all categories are included
def calculate_proportions(column):
    # Create an array of zeros for all possible nucleotides (1-4)
    props = np.zeros(4)
    # Count occurrences
    counts = np.bincount(column)[1:5] if max(column) < 5 else np.bincount(column, minlength=5)[1:5]
    # Fill available counts
    props[:len(counts)] = counts
    # Normalize
    return props / np.sum(props) if np.sum(props) > 0 else props


# Function to find dominant values
def find_dominant_values(column, threshold):
    dominant_values = np.where(column > threshold)[0]
    return dominant_values if dominant_values.size > 0 else np.nan

# Read MOI data for transfer size estimation
MOI_results = pd.read_csv(wd + "results/tables/moi_cleaned.tsv", sep="\t")
base_transfer_sizes = dict(zip(MOI_results[MOI_results['exp_line'] == exp_line]['passage'], MOI_results[MOI_results['exp_line'] == exp_line]['btk_size']))
base_transfer_sizes[0] = base_transfer_sizes[list(base_transfer_sizes.keys())[0]]

# Set simulation parameters
print(exp_line, flush = True)
genome_length, initial_population = 2100, 400
mutation_rate = 2e-5
seq_sampling_frac = 100
# total_generations = int(generation_time)
sample_nr = int(sample_nr)
# print(total_generations, flush = True)
# bottleneck_intervals, sampling_freq = 2, 2
bottleneck_intervals = int(bottleneck_freq)
if bottleneck_intervals == 2:
    R0 = 44
elif bottleneck_intervals == 3:
    R0 = 12
seq_sampling_freq = int(seq_sampling_freq)

# Initialize population as an integer NumPy array
np.random.seed(2+sample_nr)

if sample_nr == 1:
    init_population = np.tile(np.random.choice([1, 2, 3, 4], genome_length, replace=True).astype(np.uint8), (initial_population, 1))
    np.save(f"{wd}results/tables/misc/neutral-seq-sim/sequences/init_population.npy", init_population)
else:
    prev_sample_nr = sample_nr-1
    init_population = np.load(f"{wd}results/tables/misc/neutral-seq-sim/sequences/init_population.npy")
    starting_population = np.load(f"{wd}results/tables/misc/neutral-seq-sim/populations/{prev_sample_nr}.npy")


# Store total mutations over time
# mutation_counts = np.zeros((total_generations // sampling_freq)-1)

# Simulation loop
if sample_nr == 1:
    population = init_population.copy()
else:
    population = starting_population.copy()


for gen in range(((sample_nr-1)*seq_sampling_freq)+1, sample_nr*seq_sampling_freq + 1):
    print(gen, flush = True)

    bottleneck_size = calculate_transfer_size(max(1, (gen + 1) // bottleneck_intervals), base_transfer_sizes)
    print(bottleneck_size, flush = True)

    # Step 1: Mutation - Apply mutation rate
    mutation_mask = np.random.rand(*population.shape) < mutation_rate
    mutation_indices = np.where(mutation_mask)
    if mutation_indices[0].size > 0:
        population[mutation_indices] = np.array([modify_base(x) for x in population[mutation_indices]])

    # Step 2: Replication - Each genome produces R0 offspring
    if gen % bottleneck_intervals == 0:
        population = population[np.repeat(np.random.choice(population.shape[0], population.shape[0], replace=True), repeats = R0)]
    else :
        population = population[np.random.choice(population.shape[0], population.shape[0]*R0, replace=True)]


    # Step 3: Sequencing
    if gen % seq_sampling_freq == 0:
        # psg = sample_nr

        # Determine variant frequency
        #proportions = np.apply_along_axis(calculate_proportions, axis=0, arr=population).astype(np.float32)
        sampled_population = population[np.random.choice(population.shape[0], round(population.shape[0]/seq_sampling_frac))]

        # Apply function to each column
        # dominant_values_per_column = np.apply_along_axis(find_dominant_values, axis=0, arr=proportions, threshold=min_threshold)+1
        # dominant_values_per_column = dominant_values_per_column.flatten()

        # Step 4: Track mutation accumulation
        # mutation_counts[psg - 1] = np.count_nonzero(dominant_values_per_column != init_population[0, :])
        # mutation_counts[psg - 1] = np.count_nonzero((dominant_values_per_column != init_population[0, :]) & (dominant_values_per_column>0))

        np.save(f"{wd}results/tables/misc/neutral-seq-sim/sequences/{sample_nr}.npy", sampled_population)
        
        # if psg > 1 and mutation_counts[psg - 1] > mutation_counts[psg - 2]:
        #     print("######", flush = True)
        #     print(np.where(dominant_values_per_column != init_population[0, :] & (dominant_values_per_column>0)), flush = True)
    
    # Step 5: Apply bottleneck every 2 generations
    if gen % bottleneck_intervals == 0:
        population = population[np.random.choice(population.shape[0], round(bottleneck_size/3), replace=False)]
        if gen % seq_sampling_freq == 0:
            np.save(f"{wd}results/tables/misc/neutral-seq-sim/populations/{sample_nr}.npy", population)
    else:
        pass




"""
# Save results
mutation_counts_df = pd.DataFrame(mutation_counts)
mutation_counts_df.to_csv(
    f"{wd}results/tables/misc/neutral-seq-sim/mut_counts/neutral_simulation_{exp_line}_{category}_{genome_length}_{total_generations}_{seq_sampling_frac}.tsv",
    sep="\t", index=False, header=False
)

proportions_df = pd.DataFrame(proportions)
proportions_df.to_csv(
    f"{wd}results/tables/misc/neutral-seq-sim/proportions/neutral_simulation_{exp_line}_{genome_length}_{total_generations}_{seq_sampling_frac}.tsv",
    sep="\t", index=False, header=False
)
"""
