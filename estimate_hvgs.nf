nextflow.enable.dsl=2

params.inputSeuratObject = ''
params.filterEssential = ''
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
    path essential_genes_path

    output:
    path 'gene_variance.csv'

    publishDir params.outputDir, mode: 'copy'

    script:
    """
    estimate_gene_variance.py ${pseudobulk_path} ${essential_genes_path}
    """
}

workflow {
    input_seurat = Channel.fromPath(params.inputSeuratObject)
    essential_genes_path = Channel.fromPath(params.filterEssential)
    pseudobulk = PSEUDOBULK(input_seurat)
    ESTIMATE_GENE_VARIANCE(pseudobulk, essential_genes_path)
}
