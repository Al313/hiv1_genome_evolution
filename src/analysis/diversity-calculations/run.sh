#!/bin/bash


timestamp=$(date +%F_%T)

job_dir="/home/amovas/scratch/.slurm/jobs/${timestamp}"
out_dir="/home/amovas/scratch/.slurm/outs/${timestamp}"

job_script="${job_dir}/${timestamp}.job"

# Create the output directory if it doesn't exist
mkdir -p ${out_dir}
mkdir -p ${job_dir}


echo "#!/bin/bash

#SBATCH -J diversity_calc
#SBATCH --array=13-20
#SBATCH -c 1
#SBATCH --time=3-00:00:00
#SBATCH --mem=15G	
#SBATCH --output=${out_dir}/%x-%j-%a.out

echo \$SLURM_ARRAY_TASK_ID

source activate ha_proj

Rscript /home/amovas/data/genome-evo-proj/src/analysis/diversity-calculations/get_diversity.R \$SLURM_ARRAY_TASK_ID

" > ${job_script}

sbatch ${job_script}



