process FUNCOTATOR {
    tag "${sample_ID}"
    label "low"

    publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true') 

    input:
    path genome
    path index
    path dict
    path annotation_resources
    tuple val(sample_ID), path(vcf), path(vcf_index)

    output:
    path("${sample_ID}.funcotated.maf")
    
    script:
    """
    gatk Funcotator \
    -V ${vcf} \
    -R ${genome} \
    --ref-version hg38 \
    --data-sources-path ${annotation_resources} \
    -O ${sample_ID}.funcotated.maf \
    -output-file-format MAF
    """
}

// annotation-default
// https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial
// seems like max heap memory is not enough