#!/bin/bash


timestamp=$(date +%F_%T)


job_dir="/home/amovas/scratch/.slurm/jobs/${timestamp}"


if [ ! -d ${job_dir} ]; then
  mkdir -p ${job_dir};
fi

job_file=${job_dir}/neutral_seq_sim_main.job


output_dir="/home/amovas/scratch/.slurm/outs/${timestamp}"

if [ ! -d ${output_dir} ]; then
  mkdir -p ${output_dir};
fi


cat > "${job_file}" <<EOF
#!/bin/bash


#SBATCH --cpus-per-task=4
#SBATCH --time=3-24:00:00
#SBATCH --mem="50GB"
#SBATCH --output=${output_dir}/neutral_seq_sim.out

# Load R module (assuming the cluster has R installed as a module)
source activate ha_proj


# parameters

bottleneck_freq=3
tot_gen_nr=\$(( 500 * bottleneck_freq ))
seq_sampling_freq=\$(( 10 * bottleneck_freq ))
tot_seq=\$(( tot_gen_nr / seq_sampling_freq ))
exp_line="MT-2_1"
tot_run=10



# job specifications
timestamp=\$(date +%F_%T)

job_dir="/home/amovas/scratch/.slurm/jobs/\${timestamp}"


if [ ! -d \${job_dir} ]; then
  mkdir -p \${job_dir};
fi

job_file=\${job_dir}/neutral_seq_sim.job




output_dir="/home/amovas/scratch/.slurm/outs/\${timestamp}"

if [ ! -d \${output_dir} ]; then
  mkdir -p \${output_dir};
fi



rm /home/amovas/data/genome-evo-proj/results/tables/neutral-seq-sim/populations/\${bottleneck_freq}/\${exp_line}/*
rm /home/amovas/data/genome-evo-proj/results/tables/neutral-seq-sim/sequences/\${bottleneck_freq}/\${exp_line}/*

# Loop over runs

for run_nr in \$(seq 1 \${tot_run}); do

# Loop over samples

for sample_nr in \$(seq 1 \${tot_seq}); do


job_file="\${job_dir}/neutral_seq_sim_\${exp_line}_\${sample_nr}.job"

mem=\$(( sample_nr + 10 ))




echo "#!/bin/bash


#SBATCH --cpus-per-task=4
#SBATCH --time=3-24:00:00
#SBATCH --mem="\${mem}GB"
#SBATCH --output=\${output_dir}/neutral_seq_sim_\${exp_line}_\${sample_nr}.out

# Load R module (assuming the cluster has R installed as a module)
source activate ha_proj

# Run the python script
stdbuf -oL ./neutral_sequence_simulation.py \${exp_line} \${tot_gen_nr} \${run_nr} \${sample_nr} \${bottleneck_freq} \${seq_sampling_freq}

" > \${job_file}

# Submit job and capture the job ID
job_submission_output=\$(sbatch \${job_file})
echo \${job_submission_output}
job_id=\$(echo \${job_submission_output} | awk '{print \$4}')
touch \${output_dir}/\${job_id}

# Path to the file youre waiting for
FILE="/home/amovas/data/genome-evo-proj/results/tables/neutral-seq-sim/populations/\${bottleneck_freq}/\${exp_line}/run\${run_nr}/\${sample_nr}.npy"

while [ ! -f "\$FILE" ]; do
echo "Waiting for \$FILE to be created..."
sleep 10  # Check every 10 seconds
done


done
done
EOF

sbatch ${job_file}

#done
