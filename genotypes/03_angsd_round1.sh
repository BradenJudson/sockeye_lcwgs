#!/bin/bash

#SBATCH --job-name=sock_angsd
#SBATCH --account=grdi_genarcc
#SBATCH --open-mode=append
#SBATCH --partition=standard
#SBATCH --output=logs/angsd_chr.%A_%a.out
#SBATCH --error=logs/angsd_chr.%A_%a.err
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --mem=110GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-161

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS
conda activate ../lcwgs_env

REGION=`cat ../sockeye_10mb_genome.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Sockeye_Salmon
GENOME="GCF_006149115.2_Oner_1.1_genomic.fna"
BAMFILES="bam_list_cov1.txt"
OUTFOLDER="04a_genotype_likelihoods_all"

#ls -1 03_alignments_dedup_clip/*.bam > 01_info_files/"$BAMFILES"

angsd -b 01_info_files/"$BAMFILES" -ref "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/sockeye_angsd_"$REGION" \
        -nThreads 8 \
        -r "$REGION" \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -SNP_pval 1e-10 \
        -minMapQ 20 \
        -minQ 20 \
        -setMinDepth 379 \
        -setMaxDepth 6000 \
        -minInd 266 -GL 1 \
        -doGlf 2 -doMaf 1 \
        -doMajorMinor 1 \
        -doCounts 1 \
        -rmTriallelic 1e-6 \
        -dumpCounts 2 \
        -minMaf 0.01 \
        -only_proper_pairs 1 \
        -doSaf 1