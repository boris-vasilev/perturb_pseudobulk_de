nextflow.enable.dsl=2

params.inputH5AD = ''
params.outputDir = '' // outputDir for gene_variance.csv

process ESTIMATE_BASE_EXPRESSION {
    time '5m'
    cpus 20
    label 'short'

    input:
    path h5ad_path

    output:
    path 'base_expression.csv'

    publishDir params.outputDir, mode: 'copy'

    script:
    """
    base_expression.py ${h5ad_path}
    """
}

workflow {
    input_h5ad = Channel.fromPath(params.inputH5AD)
    ESTIMATE_BASE_EXPRESSION(input_h5ad)
}
