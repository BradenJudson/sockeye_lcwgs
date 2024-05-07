#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/sockeye_lcWGS
conda activate ../lcwgs_env

export PATH="$HOME/gpfs/fs7/grdi/genarcc/wp3/judsonb/software/pcangsd/bin/:$PATH"
HEADFILE=`ls -1 *.beagle.gz | tail -n 1`

cd 04a_genotype_likelihoods_all
zcat "$HEADFILE" | head -n 1 > header_beagle.txt
cat header_beagle.txt <(gunzip -c *beagle.gz | grep -v "marker") | gzip > angsd_c.beagle.gz

pcangsd -b angsd_c.beagle.gz -o ../pca_full -t 64










