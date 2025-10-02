nextflow.enable.dsl=2

params.inputSeuratObject = ''
params.outputDir = '' // outputDir for gene_variance.csv

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

process ESTIMATE_GENE_VARIANCE {
    time '10m'
    cpus 20
    label 'short'

    input:
    path pseudobulk_path

    output:
    path 'gene_variance.csv'

    publishDir params.outputDir, mode: 'copy'

    script:
    """
    estimate_gene_variance.py ${pseudobulk_path}
    """
}

workflow {
    input_seurat = Channel.fromPath(params.inputSeuratObject)
    pseudobulk = PSEUDOBULK(input_seurat)

}
