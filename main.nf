nextflow.enable.dsl=2

params.inputSeuratObject = '/rds/user/biv22/hpc-work/AMD/data/perturb/seurat/GSE264667_jurkat_normalized_singlecell_01.h5seurat'
params.numWorkers = 16
params.outputDir = 'DEGs' // outputDir for DEGs

process PSEUDOBULK {
    time '20m'
    cpus 20
    label 'short'

    input:
    path h5ad_path

    output:
    path 'pseudobulk.h5ad'

    script:
    """
    pseudobulk.py ${h5ad_path}
    """
}

process SPLIT_SEURAT_OBJECT {
    time '20m'
    cpus 20
    label 'short'

    input:
    path pseudobulk_path

    output:
    path 'pseudobulk_chunk_*.h5ad'

    script:
    """
    chunk.py ${pseudobulk_path} ${params.numWorkers}
    """
}

process DIFFERENTIAL_EXPRESSION {
    time '8h'
    cpus 2
    executor 'slurm'
    array params.numWorkers

    input:
    path seurat_chunk

    output:
    path '*.tsv', optional: true

    publishDir params.outputDir, mode: 'copy'

    script:
    """
    differential_expression_chunks.R ${seurat_chunk} ${params.outputDir}
    """
}

workflow {
    input_seurat = Channel.fromPath(params.inputSeuratObject)
    pseudobulk = PSEEUDOBULK(input_seurat)
    chunks = SPLIT_SEURAT_OBJECT(pseudobulk)
    DIFFERENTIAL_EXPRESSION(chunks.flatten())
}
