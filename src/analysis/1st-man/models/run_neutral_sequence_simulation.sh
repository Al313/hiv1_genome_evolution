#!/bin/bash



timestamp=$(date +%F_%T)

job_dir="/home/amovas/scratch/.slurm/jobs/${timestamp}"


if [ ! -d ${job_dir} ]; then
  mkdir -p ${job_dir};
fi

job_file=${job_dir}/neutral_seq_sim.job




output_dir="/home/amovas/scratch/.slurm/outs/${timestamp}"

if [ ! -d ${output_dir} ]; then
  mkdir -p ${output_dir};
fi

# Define an array of cell lines
exp_lines=("MT-2_1" "MT-2_2" "MT-4_1" "MT-4_2")


# Check if the mutation type argument is provided
if [ $# -eq 0 ]; then
  echo "Error: Mutation type is required. Please provide a mutation type (e.g., Majority or Fixed)."
  exit 1
fi


# Define an array of mutations
# mutation_types=("Majority" "Fixed")

# Parse command-line arguments for mutation types
mutation_types=("$@")


# Loop over each cell line and submit a separate job
for exp_line in "${exp_lines[@]}"; do
for mut_cat in "${mutation_types[@]}"; do


# Set the generation time based on the mutation type
if [ "$mut_cat" == "Majority" ]; then
      generation_time=180
      mem="200G"
elif [ "$mut_cat" == "Fixed" ]; then
      generation_time=1000
      mem="300G"
fi

job_file="${job_dir}/neutral_seq_sim_${exp_line}_${mut_cat}.job"



echo "#!/bin/bash


#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=${mem}
#SBATCH --output=${output_dir}/neutral_seq_sim_${exp_line}_${mut_cat}.out

# Load R module (assuming the cluster has R installed as a module)
source activate ha_proj

# Run the R script
Rscript ./neutral_sequence_simulation.R ${exp_line} ${mut_cat} ${generation_time}

" > ${job_file}

sbatch ${job_file}

done
done

