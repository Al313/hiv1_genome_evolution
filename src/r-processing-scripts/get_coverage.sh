#!/bin/bash


wd="/home/ali313/Desktop/ltee_raw/test-snakemake/freezed-data"

exps=("III" "III" "III" "III" "IV" "IV" "IV" "IV")
lines=(13 14 15 16 17 18 19 20)
exp_iii_lines=(13 14 15 16)
exp_iv_lines=(17 18 19 20)


start_pos=(500 2000 3500 5000 6500 8000)  
end_pos=(1999 3499 4999 6499 7999 9499)


for i in {0..0}; do #0..7

my_exp=${exps[${i}]}
my_line=${lines[${i}]}

bam_dir="${wd}/HIV_LTE_mappings/${my_exp}/${my_line}"

echo ${bam_dir}

for k in {0..5};do

echo ${k}

start=${start_pos[${k}]}
end=${end_pos[${k}]}



for item in ${exp_iii_lines[@]} ; do
    [[ ${my_line} == "$item" ]] && cut_pos=14
done


for item in ${exp_iv_lines[@]} ; do
    [[ ${my_line} == "$item" ]] && cut_pos=13
done

echo $cut_pos

for j in ${bam_dir}/*.bam; do echo "$(basename ${j})" | awk "{print substr(\$0,${cut_pos},3)}" && samtools depth -a -r AF324493.2:${start}-${end} ${j}  |  awk '{sum+=$3} END { print sum/NR}'; done > ${wd}/covergae/${my_line}_${k}.txt

done

done






