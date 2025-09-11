#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 8:00:00

#. /etc/profile.d/modules.sh

#module load miniconda/3
source ~/.bashrc
conda activate env_nf

nextflow main.nf \
    --inputSeuratObject /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/AMD/data/perturb/seurat/rpe1_normalized_singlecell_01.h5seurat \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/AMD/data/perturb/DEGs/RPE1 \
    --numWorkers 50
