#!/bin/bash

. /etc/profile.d/modules.sh
module load gcc/11

#module load miniconda/3
source ~/.bashrc

conda activate env_nf

echo "Estimating K562_essential Base Expression"
nextflow estimate_base_expression.nf -resume \
   --inputH5AD /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/K562_essential_raw_singlecell_01.h5ad \
   --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/K562_essential

# echo "Estimating K562_GenomeWide Base Expression"
# nextflow estimate_base_expression.nf -resume \
#    --inputH5AD /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/K562_gwps_raw_singlecell.h5ad \
#    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/K562_GenomeWide \
#    --filterEssential /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/essential_genes.txt

echo "Estimating Jurkat Base Expression"
nextflow estimate_base_expression.nf -resume \
    --inputH5AD /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/jurkat_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/Jurkat

echo "Estimating HepG2 Base Expression"
nextflow estimate_base_expression.nf -resume \
    --inputH5AD /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/hepg2_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/HepG2

echo "Estimating RPE1 Base Expression"
nextflow estimate_base_expression.nf -resume \
    --inputH5AD /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/rpe1_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/RPE1
