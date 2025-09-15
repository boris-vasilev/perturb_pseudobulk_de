nextflow.enable.dsl=2

params.inputSeuratObject = '/rds/user/biv22/hpc-work/AMD/data/perturb/seurat/GSE264667_jurkat_normalized_singlecell_01.h5seurat'
params.numWorkers = 16
params.outputDir = 'DEGs' // outputDir for DEGs

process SPLIT_SEURAT_OBJECT {
    memory '50GB'
    time '6h'

    input:
    path h5ad_path

    output:
    path 'pseudobulk_chunk_*.h5ad'

    // pseudobulk only 1 chunk (full data)
    script:
    """
    pseudobulk.py $h5ad_path 1
    """
}

process CHUNKS_DIFFERENTIAL_EXPRESSION {
    memory '50GB'
    time '8h'
    cpus 16
    //maxForks 8

    input:
    path seurat_chunk

    output:
    path '*.tsv', optional: true

    publishDir params.outputDir, mode: 'copy'

    script:
    """
    DE_random_effect_batch.R ${seurat_chunk} ${params.outputDir}
    """
}

workflow {
    input_seurat = Channel.fromPath(params.inputSeuratObject)
    input_seurat.view()
    seurat_chunks = SPLIT_SEURAT_OBJECT(input_seurat)
    CHUNKS_DIFFERENTIAL_EXPRESSION(seurat_chunks.flatten())
}
