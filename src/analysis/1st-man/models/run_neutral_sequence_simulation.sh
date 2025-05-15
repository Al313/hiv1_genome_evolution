#!/bin/bash



# parameters

# Define an array of cell lines
exp_lines=("MT-2_1")

# Check if the mutation type argument is provided
if [ $# -eq 0 ]; then
  echo "Error: Mutation type is required. Please provide a mutation type (e.g., Majority or Fixed)."
  exit 1
fi

# Define an array of mutations
# mutation_types=("Majority" "Fixed")

# Parse command-line arguments for mutation types
mutation_types=("$@")


# Set the generation time based on the mutation type
if [ "$mutation_types" == "Majority" ]; then
      tot_gen_nr=100
      mem="50G"
elif [ "$mutation_types" == "Fixed" ] || [ "$mutation_types" == "Minority" ]; then
      tot_gen_nr=100
      mem="50G"
fi

bottleneck_freq=2
seq_sampling_freq=20






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




# Loop over each cell line and submit a separate job
for exp_line in "${exp_lines[@]}"; do
for mut_cat in "${mutation_types[@]}"; do

tot_seq=$(( tot_gen_nr / (seq_sampling_freq) ))

for sample_nr in $(seq 1 ${tot_seq}); do


job_file="${job_dir}/neutral_seq_sim_${exp_line}_${mut_cat}_${sample_nr}.job"



echo "#!/bin/bash


#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=${mem}
#SBATCH --output=${output_dir}/neutral_seq_sim_${exp_line}_${mut_cat}_${sample_nr}.out

# Load R module (assuming the cluster has R installed as a module)
source activate ha_proj

# Run the python script
stdbuf -oL ./neutral_sequence_simulation.py ${exp_line} ${mut_cat} ${tot_gen_nr} ${sample_nr} ${bottleneck_freq} ${seq_sampling_freq}

" > ${job_file}

# Submit job and capture the job ID
job_submission_output=$(sbatch ${job_file})
echo ${job_submission_output}
job_id=$(echo ${job_submission_output} | awk '{print $4}')
touch ${output_dir}/${job_id}

# Path to the file you're waiting for
FILE="/home/amovas/data/genome-evo-proj/results/tables/misc/neutral-seq-sim/populations/${sample_nr}.npy"

while [ ! -f "$FILE" ]; do
echo "Waiting for $FILE to be created..."
sleep 10  # Check every 10 seconds
done


done
done
done
