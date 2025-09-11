#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 4:00:00

#. /etc/profile.d/modules.sh

#module load miniconda/3
. /etc/profile.d/modules.sh
module load gcc/11

source ~/.bashrc
conda activate env_nf

nextflow main.nf \
    --inputSeuratObject /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/AMD/data/perturb/seurat/K562_essential_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/AMD/data/perturb/DEGs/K562_pseudo \
    --numWorkers 75 

