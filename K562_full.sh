#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 24:00:00

#. /etc/profile.d/modules.sh

#module load miniconda/3
source ~/.bashrc
conda activate env_nf

nextflow full.nf \
    --inputSeuratObject /rds/user/biv22/hpc-work/AMD/data/perturb/seurat/K562_gwps_normalized_singlecell_01.h5seurat \
    --outputDir /rds/user/biv22/hpc-work/AMD/data/perturb/DEGs/K562
