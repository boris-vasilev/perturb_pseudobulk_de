#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 10:00:00

#. /etc/profile.d/modules.sh

#module load miniconda/3
source ~/.bashrc
conda activate env_nf

nextflow main.nf \
    --inputSeuratObject /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/AMD/data/perturb/seurat/GSE264667_jurkat_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/AMD/data/perturb/DEGs/Jurkat_pseudo \
    --numWorkers 75
