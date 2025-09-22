#!/bin/bash

. /etc/profile.d/modules.sh
module load gcc/11

#module load miniconda/3
source ~/.bashrc

nextflow main.nf -resume \
    --inputSeuratObject /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/rpe1_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb/RPE1 \
    --numWorkers 100
