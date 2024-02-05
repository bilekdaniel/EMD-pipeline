Channel
    .fromPath(params.sampleCSV)
    .splitCsv(header: true)
    .map { row -> [sample_ID: row.names, fq1: "${params.samples}${row.fq1}", fq2: "${params.samples}${row.fq2}"] }
    .set { inv_samples }

process DESEQ2 {
    publishDir(path: "${results}/deseq2_ndmm_mrd", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), path(quant_folder)

    output: 
    path("results")

    """
    Rscript /scratch/project/open-27-18/mrd_trans/scripts/deseq2_contrast.R ${params.sampleCSV} ${params.deseq_run_name} ${baseDir}
    """
}

