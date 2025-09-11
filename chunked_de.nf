nextflow.enable.dsl=2

params.chunksDir = "/rds/user/biv22/hpc-work/AMD/data/perturb/seurat_Jurkat_chunks/seurat_chunk_*.RDS"
params.numWorkers = 16
params.outputDir = 'DEGs' // outputDir for DEGs

process CHUNKS_DIFFERENTIAL_EXPRESSION {
    memory '10GB'
    time '2h'
    cpus 1

    input:
    path seurat_chunk

    output:
    path '*.tsv'

    publishDir params.outputDir, mode: 'copy'

    script:
    """
    differential_expression_chunks.R ${seurat_chunk}
    """
}

workflow {
    seurat_chunks = Channel.fromPath(params.chunksDir)
    seurat_chunks.view { "File found: $it" } // This will print each file found 
    CHUNKS_DIFFERENTIAL_EXPRESSION(seurat_chunks)
}


