process ADD_READ_GROUP {
    label: "low"
    publishDir(path: "${params.output}${sample_ID}/star", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), path(BAM_FILE)

    output:
    tuple val(sample_ID), path ("${sample_ID}_rgAdded_aligned.sortedByCoord.out.bam"), emit: rgBAM_FILE

    """
    picard AddOrReplaceReadGroups --RGLB lib1 --RGPL illumina --RGPU ${sample_ID} -SM ${sample_ID} \
    --VALIDATION_STRINGENCY LENIENT -I ${BAM_FILE} -O ${sample_ID}_rgAdded_aligned.sortedByCoord.out.bam 
    """
}