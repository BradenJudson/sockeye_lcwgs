#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS
conda activate ../lcwgs_env

export PATH="$HOME/gpfs/fs7/grdi/genarcc/wp3/judsonb/software/pcangsd/bin/:$PATH"

cd 04b_genotype_likelihoods_200kSNPs
zcat sockeye_angsd200k_NC_042541.1:50000000-58386872.beagle.gz | head -n 1 > header_beagle.txt
cat header_beagle.txt <(gunzip -c *beagle.gz | grep -v "marker") | gzip > angsd_c.beagle.gz

cd ..
pcangsd -b 04b_genotype_likelihoods_200kSNPs/angsd_c.beagle.gz -o pca_200k -t 64


