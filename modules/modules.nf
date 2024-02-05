
process MULTIQC { 
    publishDir(path: "${params.output}MultiQC", mode: 'copy', overwrite: 'true')
    
    input:
    tuple val(sample_ID), path(star_qualimap)
    tuple val(sample_ID), path(read_1_no_rrna)
    tuple val(sample_ID), path(read_2_no_rrna)
    tuple val(sample_ID), path(quant_folder)

    output:
    path "./multiqc_report.html"

    """
    multiqc .
    """

}
process FASTQC_BOTH {

    publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), path(read_1_no_rrna)
    tuple val(sample_ID), path(read_2_no_rrna)

    output:
    path "${sample_ID}_1_ribo_norrna_fastqc.html", optional: true
    path "${sample_ID}_1_ribo_norrna_fastqc.zip", optional: true
    path "${sample_ID}_2_ribo_norrna_fastqc.html", optional: true
    path "${sample_ID}_2_ribo_norrna_fastqc.zip", optional: true

    """
    fastqc -t 10 ${read_1_no_rrna} ${read_2_no_rrna}
    """

}
// process trinity {

//     publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true')

//     input:
//     tuple val(sample_ID), path(read_1_no_rrna) //perhaps has to be .gz
//     tuple val(sample_ID), path(read_2_no_rrna)


//     output:
//     tuple val(sample_ID), path("./trinity_out_dir"), optional: true

//     """
//     Trinity --seqType fq \
//     --left ${read_1_no_rrna} --right ${read_2_no_rrna} \
//     --max_memory 250G \
//     --CPU 32 \
//     """
// }



