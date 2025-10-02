#!/bin/bash

. /etc/profile.d/modules.sh
module load gcc/11

#module load miniconda/3
source ~/.bashrc

echo "Estimating K562_essential HGVs"
nextflow estimate_hgvs.nf \
    --inputSeuratObject /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/K562_essential_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/K562_essential

echo "Estimating K562_GenomeWide HGVs"
nextflow estimate_hgvs.nf \
    --inputSeuratObject /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/K562_gwps_raw_singlecell.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/K562_GenomeWide

echo "Estimating Jurkat HGVs"
nextflow estimate_hgvs.nf \
    --inputSeuratObject /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/jurkat_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/Jurkat

echo "Estimating HepG2 HGVs"
nextflow estimate_hgvs.nf \
    --inputSeuratObject /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/hepg2_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/HepG2

echo "Estimating RPE1 HGVs"
nextflow estimate_hgvs.nf \
    --inputSeuratObject /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_pseudobulk_de/data/rpe1_raw_singlecell_01.h5ad \
    --outputDir /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y/biv22/data/perturb_QC/RPE1