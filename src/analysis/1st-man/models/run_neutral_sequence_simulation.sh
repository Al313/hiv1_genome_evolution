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



echo "#!/bin/bash


#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=30G
#SBATCH --output=${output_dir}/neutral_seq_sim.out

# Load R module (assuming the cluster has R installed as a module)
source activate ha_proj

# Run the R script
Rscript ./neutral_sequence_simulation.R

" > ${job_file}

sbatch ${job_file}




