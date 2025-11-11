#!/usr/bin/env python3




import numpy as np
import pandas as pd
import sys
import os
import math
import re

# Get command-line arguments
args = sys.argv[1:] if len(sys.argv) > 1 else ["MT-2_1", "Majority", 180]
exp_line, generation_time, sample_nr, bottleneck_freq, seq_sampling_freq = args

# Determine working directory
if os.path.exists("/home/amovas/"):
    print("Remote HPC Connection!")
    wd = "/home/amovas/data/genome-evo-proj/"
else:
    print("Local PC Connection!")
    wd = "/Users/alimos313/Documents/studies/phd/hpc-research/genome-evo-proj/"
    exp_line, generation_time = "MT-2_1", 180



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


# Read MOI data for transfer size estimation
MOI_results = pd.read_csv(wd + "results/tables/moi_cleaned.tsv", sep="\t")
base_transfer_sizes = dict(zip(MOI_results[MOI_results['exp_line'] == exp_line]['passage'], MOI_results[MOI_results['exp_line'] == exp_line]['btk_size']))
base_transfer_sizes[0] = base_transfer_sizes[list(base_transfer_sizes.keys())[0]]

# Set simulation parameters
print(exp_line, flush = True)

if re.match(r"^MT-2(?:_|$)", exp_line):
    neutral_pos = 2000
elif re.match(r"^MT-4(?:_|$)", exp_line):
    neutral_pos = 1500

print(exp_line)
print(neutral_pos)

initial_population = 400
mutation_rate = 2e-5
seq_sampling_frac = 1000
infection_success_rate = 0.1
sample_nr = int(sample_nr)
bottleneck_intervals = int(bottleneck_freq)
if bottleneck_intervals == 2:
    R0 = 44
elif bottleneck_intervals == 3:
    R0 = 12
seq_sampling_freq = int(seq_sampling_freq)




# Define path and filename
folder_path_population = f"{wd}results/tables/neutral-seq-sim/populations/{bottleneck_intervals}/{exp_line}"
folder_path_sequence = f"{wd}results/tables/neutral-seq-sim/sequences/{bottleneck_intervals}/{exp_line}"
# Ensure the directory exists
os.makedirs(folder_path_population, exist_ok=True)
os.makedirs(folder_path_sequence, exist_ok=True)


# Create and save initial ancestral population and read previous populations
np.random.seed(sample_nr+2)

if sample_nr == 1:
    init_population = np.tile(np.random.choice([1, 2, 3, 4], neutral_pos, replace=True).astype(np.uint8), (initial_population, 1))
    np.save(f"{folder_path_sequence}/init_population.npy", init_population)
else:
    prev_sample_nr = sample_nr-1
    init_population = np.load(f"{folder_path_sequence}/init_population.npy")
    starting_population = np.load(f"{folder_path_population}/{prev_sample_nr}.npy")






# Simulation loop
if sample_nr == 1:
    population = init_population.copy()
else:
    population = starting_population.copy()


for gen in range(((sample_nr-1)*seq_sampling_freq)+1, sample_nr*seq_sampling_freq + 1):
    print(gen, flush = True)

    bottleneck_size = calculate_transfer_size(max(1, (gen + 1) // bottleneck_intervals), base_transfer_sizes)
    bottleneck_size = round(bottleneck_size*infection_success_rate)
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

        # Sampling for sequencing
        sampled_population = population[np.random.choice(population.shape[0], round(population.shape[0]/seq_sampling_frac))]


        np.save(f"{folder_path_sequence}/{sample_nr}.npy", sampled_population)
        
    
    # Step 5: Apply bottleneck every 2 generations
    if gen % bottleneck_intervals == 0:
        population = population[np.random.choice(population.shape[0], bottleneck_size, replace=False)]
        if gen % seq_sampling_freq == 0:
            np.save(f"{folder_path_population}/{sample_nr}.npy", population)
    else:
        pass




