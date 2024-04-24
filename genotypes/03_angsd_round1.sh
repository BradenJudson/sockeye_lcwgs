#!/bin/bash

#SBATCH --job-name=sock_angsd
#SBATCH --account=grdi_genarcc
#SBATCH --open-mode=append
#SBATCH --partition=standard
#SBATCH --output=logs/angsd_chr.%A_%a.out
#SBATCH --error=logs/angsd_chr.%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-161

source ~/.bashrc

REGION=`cat ../sockeye_genome_10mb_regions.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Sockeye_Salmon
GENOME="GCF_006149115.2_Oner_1.1_genomic.fna"
BAMFILES="bam_list.txt"
OUTFOLDER="04_genotype_likelihoods"

ls -1 03_alignments_dedup_clip/*.bam > 01_info_files/"$BAMFILES"

angsd -b 01_info_files/"$BAMFILES" -ref "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/sockeye_angsd_"$REGION" \
        -nThreads 8 \
        -r "$REGION" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -SNP_pval 1e-6  \
        -minMapQ 20 \
        -minQ 20 \
        -setMinDepth 300 \
        -setMaxDepth 6000 \
        -minInd 300 -GL 1 \
        -doMaf 1 -doMajorMinor 1 \
        -doCounts 1 -skipTriallelic 1