#!/bin/bash




timestamp=$(date +%F_%T)

job_dir="/home/amovas/scratch/.slurm/jobs/${timestamp}"


if [ ! -d ${job_dir} ]; then
  mkdir -p ${job_dir};
fi

job_file=${job_dir}/snakemake.job




output_dir="/home/amovas/scratch/.slurm/outs/${timestamp}"

if [ ! -d ${output_dir} ]; then
  mkdir -p ${output_dir};
fi



echo "#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --output=${output_dir}/snakemake.out

source activate ha_proj

snakemake \
        -s mysnakefile.smk \
	--jobs 456 \
        --default-resource mem_mb=50000 \
        --cluster '
        sbatch \
                --cpus-per-task {threads} \
                --mem {resources.mem_mb} \
                --time 2:0:0 \
                --output logs/snake.out \
                --error logs/snake.err'
" > ${job_file}

sbatch ${job_file}




