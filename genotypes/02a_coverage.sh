#!/bin/bash

#SBATCH --job-name=lcwgs_cov
#SBATCH --mem-per-cpu=5G
#SBATCH --account=grdi_genarcc
#SBATCH --time=48:00:00


source ~/.bashrc

for file in *.bam;
do
        sample="$file"
        avcove=$(samtools depth -a "$file" | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR}')
        echo ${sample} ${avcove}
done > ../01_info_files/bam_coverage.txt