include { TRIM_GALORE_PE                                } from './modules/trim_galore'
include { RIBO_DET                                      } from './modules/ribo_det'
include { SALMON                                        } from './modules/DE/salmon'
include { DESEQ2                                        } from './modules/DE/deseq2'

Channel
    .fromPath(params.sampleCSV)
    .splitCsv(header: true)
    .map { row -> [sample_ID: row.names, fq1: "${params.samples}${row.fq1}", fq2: "${params.samples}${row.fq2}"] }
    .set { inv_samples }

workflow {    

    TRIM_GALORE_PE(inv_samples)

    RIBO_DET(TRIM_GALORE_PE.out)

    SALMON(RIBO_DET.out.reads_wo_rrna1, RIBO_DET.out.reads_wo_rrna2) | collect | DESEQ2

}