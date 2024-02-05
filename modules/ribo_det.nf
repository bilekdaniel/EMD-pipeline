process RIBO_DET {
    label "ribo_det"
    tag "${sample_ID}"
    publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), path(reads_1)
    tuple val(sample_ID), path(reads_2)

    output:
    tuple val(sample_ID), path("${sample_ID}_1_ribo_norrna.fastq.gz"), emit: reads_wo_rrna1
    tuple val(sample_ID), path("${sample_ID}_2_ribo_norrna.fastq.gz"), emit: reads_wo_rrna2
    path "${sample_ID}_1_ribo_rrna.fastq", optional: true
    path "${sample_ID}_2_ribo_rrna.fastq", optional: true

    """
    ribodetector_cpu -t ${task.cpus} -l 150 --chunk_size 1800 -i ${reads_1} ${reads_2} -e rrna \
    -o ${sample_ID}_1_ribo_norrna.fastq.gz ${sample_ID}_2_ribo_norrna.fastq.gz \
    -r ${sample_ID}_1_ribo_rrna.fastq ${sample_ID}_2_ribo_rrna.fastq
    """
}