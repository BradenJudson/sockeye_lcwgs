#!/bin/bash

INPUT="02_raw_reads"

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS

ls "$INPUT"/*R1.fastq.gz | perl -pe 's/_R1.*//g' | perl -pe "s/$INPUT\///g" > 01_info_files/samples.txt