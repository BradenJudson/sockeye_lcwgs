#!/bin/bash


#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --job-name=sock_popSFS
#SBATCH --output=logs_popSFS/popSFS.out
#SBATCH --error=logs_popSFS/popSFS.err
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=64


cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS
conda activate ../lcwgs_env

POPLIST=`ls -1 05_pop_bams/text_files | cut -d_ -f1`
GENOMEFOLDER=/gpfs/fs7/grdi/genarcc/common/genomes/Sockeye_Salmon
GENOME="GCF_006149115.2_Oner_1.1_genomic.fna"
OUTFOLDER="05_pop_bams/sfs_outputs"
SITES="01_info_files/sockeye_angsd_snplist.txt"

dos2unix 05_pop_bams/text_files/*.txt

for POP in $POPLIST
do
        echo $POP
        angsd -b 05_pop_bams/text_files/"$POP"_bams.txt \
        -ref "$GENOMEFOLDER"/"$GENOME" -anc "$GENOMEFOLDER"/"$GENOME" \
        -uniqueOnly 1 -out "$OUTFOLDER"/"$POP" -C 50 \
        -minMapQ 20 -minQ 20 -GL 1 -doSaf 1 -doCounts 1 \
        -only_proper_pairs -remove_bads 1 \
        -sites "$SITES"
done