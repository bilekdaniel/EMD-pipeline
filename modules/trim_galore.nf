process TRIM_GALORE_PE {
    label "trim_galore"
    tag "${sample_ID}"
    publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), val(fq1), val(fq2)

    output:
    tuple val(sample_ID), path("${sample_ID}_1_val_1.fq.gz")
    tuple val(sample_ID), path("${sample_ID}_2_val_2.fq.gz")

    """
    trim_galore --illumina -q 20 --clip_R2 3 --illumina --paired ${fq1} ${fq2} --cores ${task.cpus}
    """
}
