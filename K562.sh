#!/bin/bash
. /etc/profile.d/modules.sh
module load gcc/11

#module load miniconda/3
source ~/.bashrc
conda activate env_nf

nextflow mega.nf -resume \
    --inputSeuratObject /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/K562_gwps_raw_singlecell.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb/K562_GenomeWide \
    --numWorkers 100

