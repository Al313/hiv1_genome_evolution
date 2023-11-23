#!/bin/bash




job_dir="/home/amovas/scratch/.slurm/jobs"
out_dir="/home/amovas/scratch/.slurm/outs"


if [ ! -d ${job_dir} ]; then
  mkdir -p ${job_dir};
fi

if [ ! -d ${out_dir} ]; then
  mkdir -p ${out_dir};
fi


job_file=${job_dir}/snakemake.job

echo "#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --output=${out_dir}/snakemake.out

source activate ha_proj

snakemake \
        -s mysnakefile.smk \
	-R extract_annotation \
        --jobs 432 \
        --default-resource mem_mb=2000 \
        --cluster '
        sbatch \
                --cpus-per-task {threads} \
                --mem {resources.mem_mb} \
                --time 2:0:0 \
                --output logs/snake.out \
                --error logs/snake.err'
" > ${job_file}

sbatch ${job_file}




