#!/bin/bash

source ~/.bashrc
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS/04a_genotype_likelihoods_all
conda activate ../../lcwgs_env

# Unzip and concatenate minor allele frequency output files.
gunzip -c *.mafs.gz | cut -f1-4 | grep -v "chromo" > ../01_info_files/angsd_snplist_2sort.txt

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS/01_info_files

# Sort SNPs by position. Retains only polymorphic sites.
sort -k1,1 -k2,2n angsd_snplist_2sort.txt > sockeye_angsd_snplist.txt

# Index corresponding snp lists.
angsd sites index sockeye_angsd_snplist.txt

# Randomly choose 200k SNPs and index.
shuf -n 200000 sockeye_angsd_snplist.txt | sort -k1,1 -k2,2n > sockeye_200k_rSNPs.txt
angsd sites index sockeye_200k_rSNPs.txt