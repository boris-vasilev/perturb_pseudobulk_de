nextflow.enable.dsl=2

params.inputSeuratObject = '/rds/user/biv22/hpc-work/AMD/data/perturb/seurat/K562_gwps_normalized_singlecell_01.h5seurat'
params.outputDir = 'DEGs' // outputDir for DEGs

process PROCESS_SEURAT_OBJECT {
    memory '500GB'
    time '24h'

    input:
    path h5seurat_path
    
    output:
    path '*.tsv'

    publishDir params.outputDir, mode: 'copy'
    
    script:
    """
    differential_expression_full.R $h5seurat_path
    """
}

workflow {
    input_seurat = Channel.fromPath(params.inputSeuratObject)
    input_seurat.view()
    PROCESS_SEURAT_OBJECT(input_seurat)
}


