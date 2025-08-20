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


#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mem=48G
#SBATCH --output=${output_dir}/ld_calc.out

# Load R module (assuming the cluster has R installed as a module)
source activate ha_proj


# parameters

passage=300


Rscript /home/amovas/data/genome-evo-proj/src/analysis/ld/get_ld.R \${passage}

EOF

sbatch ${job_file}


