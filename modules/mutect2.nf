process MUTECT2 {
  tag "${sample_ID}"
  label 'mutect2'

 input: 
    path genome
    path index
    path dict
    path POM
    path POM_INDEX
    path germline_resources
    path germline_resources_index
    tuple val(sample_ID), path(bam), path(bai)


  output:
    tuple val(sample_ID), path("${sample_ID}_unfiltered.vcf.gz"), path("${sample_ID}_unfiltered.vcf.gz.stats"), path("${sample_ID}_unfiltered.vcf.gz.tbi"), path("read-orientation-model.tar.gz")

  script: 
  """
  gatk Mutect2 \
  -R ${genome} \
  -I ${bam} \
  --germline-resource ${germline_resources} \
  --panel-of-normals ${POM} \
  --f1r2-tar-gz f1r2.tar.gz \
  -O ${sample_ID}_unfiltered.vcf.gz


  gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz
  """
}