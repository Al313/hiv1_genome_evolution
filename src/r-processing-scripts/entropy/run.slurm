#!/bin/bash


timestamp=$(date +%F_%T)

#SBATCH -J hello
#SBATCH --array=13-20
#SBATCH -c 1
#SBATCH --time=2-24:00:00
#SBATCH --mem=10G
#SBATCH --output=/home/amovas/scratch/.slurm/outs/${timestamp}/%j-%a.out

echo $SLURM_ARRAY_TASK_ID

source activate ha_proj

Rscript get_shannon.R $SLURM_ARRAY_TASK_ID
