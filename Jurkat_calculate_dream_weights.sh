#!/bin/bash
#SBATCH -A MRC-BSU-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH -N 1
#SBATCH --cpus-per-task 16
#SBATCH -t 10:00:00

. /etc/profile.d/modules.sh

module load gcc/11

#module load miniconda/3
source ~/.bashrc
conda activate env_nf


Rscript calculate_dream_weights.R
