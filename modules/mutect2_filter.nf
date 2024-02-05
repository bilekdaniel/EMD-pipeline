process MUTECT2_FILTER {
  tag "${sample_ID}"
  label 'medium'
  publishDir(path: "${params.output}${sample_ID}/MUTECT2_FILTER", mode: 'copy', overwrite: 'true')

 input: 
    path genome
    path index
    path dict
    tuple val(sample_ID), path(vcf), path(stats), path(vcf_dict), path(read_model)


  output:
    tuple val(sample_ID), path("${sample_ID}_filtered.vcf.gz"), path("${sample_ID}_filtered.vcf.gz.tbi")

  script: 
  """
    gatk FilterMutectCalls \
    -R ${genome} \
    -V ${vcf} \
    --orientation-bias-artifact-priors ${read_model} \
    -O ${sample_ID}_filtered.vcf.gz
  """
}

