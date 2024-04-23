#!/bin/bash

#SBATCH --job-name=lcwgs_dedup
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc

source ~/.bashrc


INPUT="03_alignments"
OUTPUT="04_alignments_deduplicated"
METRICS="11_metrics"
JAVA_OPTS="-Xmx80G"
TMPDIR="99_tmp"


# Remove duplicate alignments
for file in $(ls "$INPUT"/*.bam | perl -pe 's/\.bam//g')
do
        name=$(basename "$file")
        echo "Deduplicating sample: $name"

        bam clipOverlap \
                --in "$INPUT"/"$name".bam \
                --out "$OUTPUT"/"$name".clip.bam \
                --stats

        samtools index "$OUTPUT"/"$name".clip.bam

done