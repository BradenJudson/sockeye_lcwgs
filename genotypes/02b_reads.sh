#!/bin/bash

#SBATCH --job-name=lcwgs_reads
#SBATCH --mem-per-cpu=5G
#SBATCH --account=grdi_genarcc
#SBATCH --time=48:00:00


source ~/.bashrc

for file in *.bam;
do
        sample="$file"
        nreads=$(samtools view -c -F 260 "$file");
        echo ${sample} ${nreads}
done > ../01_info_files/bam_reads.txt