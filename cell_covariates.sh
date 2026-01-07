#!/bin/bash

. /etc/profile.d/modules.sh
module load gcc/11

#module load miniconda/3
source ~/.bashrc

conda activate env_nf

echo "K562_essential cell covariates"
python bin/cell_covariates.py \
    /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/K562_essential_raw_singlecell_01.h5ad \
    K562_essential

echo "K562_GenomeWide cell covariates"
python bin/cell_covariates.py \
    /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/K562_gwps_raw_singlecell.h5ad \
    K562_GenomeWide

echo "Jurkat cell covariates"
python bin/cell_covariates.py \
    /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/jurkat_raw_singlecell_01.h5ad \
    Jurkat

echo "HepG2 cell covariates"
python bin/cell_covariates.py \
    /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/hepg2_raw_singlecell_01.h5ad \
    HepG2

echo "RPE1 cell covariates"
python bin/cell_covariates.py \
    /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/rpe1_raw_singlecell_01.h5ad \
    RPE1