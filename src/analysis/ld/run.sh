#!/bin/bash


timestamp=$(date +%F_%T)


job_dir="/home/amovas/scratch/.slurm/jobs/${timestamp}"


if [ ! -d ${job_dir} ]; then
  mkdir -p ${job_dir};
fi



output_dir="/home/amovas/scratch/.slurm/outs/${timestamp}"

if [ ! -d ${output_dir} ]; then
  mkdir -p ${output_dir};
fi



exps=(iii)
exp_lines=(15)

count=0
max=200000

for exp in "${exps[@]}"; do
 for exp_line in "${exp_lines[@]}"; do
  for item in /home/amovas/shared/genome-evo-proj/data/processed-data/mappings/pipeline-outputs/${exp}/${exp_line}/15MT4EXPIIIVP390seq18062021_S18_L001_sorted.bam; do
 
  job_file=${job_dir}/${count}.job
  echo "$item"
  ((count++))

  if [ "$count" -ge "$max" ]; then
   exit 0
  fi

  cat > "${job_file}" <<EOF
#!/bin/bash


#SBATCH --cpus-per-task=30
#SBATCH --time=24:00:00
#SBATCH --mem=50G
#SBATCH --output=${output_dir}/ld_calc_${count}.out

# Load R module (assuming the cluster has R installed as a module)
source activate ha_proj


Rscript /home/amovas/data/genome-evo-proj/src/analysis/ld/get_ld.R ${item}

EOF

  sbatch ${job_file}

  done
 done
done

