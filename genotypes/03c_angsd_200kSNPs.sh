#!/bin/bash

#SBATCH --account=grdi_genarcc
#SBATCH --open-mode=append
#SBATCH --partition=standard
#SBATCH --job-name=angsd200k_sockeye
#SBATCH --output=logs_r2/angsd200k.out
#SBATCH --error=logs_r2/angsd200k.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=1

source ~/.bashrc

#REGION=`cat ../sockeye_genome_10mb_regions.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Sockeye_Salmon
GENOME="GCF_006149115.2_Oner_1.1_genomic.fna"
BAMFILES="bam_list.txt"
OUTFOLDER="04b_genotype_likelihoods_200kSNPs"
SITES="01_info_files/sockeye_200k_rSNPs.txt"

ls -1 03_alignments_dedup_clip/*.bam > 01_info_files/"$BAMFILES"

angsd sites index "$SITES"

angsd -b 01_info_files/"$BAMFILES" -ref "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/sockeye_angsd200k \
        -nThreads 8 \
        -setMinDepth 300 \
        -setMaxDepth 7000 \
        -minInd 300 \
        -GL 1 -doGlf 2 -doMaf 1 \
        -doMajorMinor 3 \
        -doCounts 1 \
        -rmTriallelic 1e-6 \
        -minMaf 0.05 \
        -minQ -minMapQ 20 \
        -sites "$SITES"