#!/bin/bash

#SBATCH --job-name=angsd_sockeye
#SBATCH --time=168:00:60
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=15G
#SBATCH --account=grdi_genarcc


source ~/.bashrc


GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Sockeye_Salmon
GENOME="GCF_006149115.2_Oner_1.1_genomic.fna"
BAMFILES="bam_list.txt"
OUTFOLDER="04_genotype_likelihoods"


ls -1 03_alignments_dedup_clip/*.bam > 01_info_files/"$BAMFILES"

angsd -b 01_info_files/"$BAMFILES" -ref "$GENOMEFOLDER"/"$GENOME" \
        -out "$OUTFOLDER"/sockeye_angsd1 \
        -nThreads 8 \
        -dobcf 1 \
        -gl 1 \
        -dopost 1 \
        -dogeno 5 \
        -doGlf 2 \
        -doMajorMinor 1 \
        -domaf 1 \
        -doCounts 1 \
        -dumpCounts 2 \
        -doQsDist 1 \
        -minMapQ 30 \
        -minQ 30 \
        -minInd 5 \
        -SNP_pval 2e-6 \
        -uniqueOnly 1 \
        -minMaf 0.05 \
        -remove_bads \
        -only_proper_pairs 1
