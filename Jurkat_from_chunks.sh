#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH -t 1:20:00

#. /etc/profile.d/modules.sh

#module load miniconda/3
source ~/.bashrc
conda activate env_nf

nextflow chunked_de.nf -resume \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/AMD/data/perturb/DEGs/Jurkat \
    --numWorkeres 16
