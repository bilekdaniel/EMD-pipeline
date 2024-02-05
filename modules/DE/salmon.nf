process SALMON {
    label: "medium"
    publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), path(reads_wo_rrna1)
    tuple val(sample_ID), path(reads_wo_rrna2)

    output:
    tuple val(sample_ID), path("salmon"), emit: quant_folder

    """
    salmon quant -i ${params.salmon_index} --gcBias \
    --reduceGCMemory --seqBias -l A -1 ${reads_wo_rrna1} -2 ${reads_wo_rrna2} \
    --validateMappings --writeUnmappedNames -p ${task.cpus} -o salmon
    """
}
