process GETPILEUP {
  tag "${sample_ID}"
  label 'mutect2'


    input:
    path genome
    path index
    path dict
    path germline_resources
    path germline_resources_index
    tuple val(sample_ID), path(bam), path(bai)

    output:
    tuple val(sample_ID), path("${sample_ID}.pileups.table")

    script:
    """
    gatk GetPileupSummaries \
    -I ${bam} \
    -R ${genome} \  
    -V ${germline_resources} \
    -L ${germline_resources} \
    -O ${sample_ID}.pileups.table
    """
}