#!/bin/bash

source ~/.bashrc
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS
conda activate ../lcwgs_env

export PATH="$HOME/gpfs/fs7/grdi/genarcc/wp3/judsonb/software/pcangsd/bin/:$PATH"

cd 04a_genotype_likelihoods_all
zcat sockeye_angsd_NC_042541.1:50000000-58386872.beagle.gz | head -n 1 > header_beagle.txt
cat header_beagle.txt <(gunzip -c *beagle.gz | grep -v "marker") | gzip > angsd_c.beagle.gz

pcangsd -b angsd_c.beagle.gz -o ../pca_full -t 64


