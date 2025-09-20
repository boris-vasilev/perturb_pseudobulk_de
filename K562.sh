#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 4:00:00

. /etc/profile.d/modules.sh
module load gcc/11

#module load miniconda/3
source ~/.bashrc
conda activate env_nf

nextflow mega.nf -resume \
    --inputSeuratObject data/K562_gwps_raw_singlecell.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb/K562_GenomeWide \
    --numWorkers 30

