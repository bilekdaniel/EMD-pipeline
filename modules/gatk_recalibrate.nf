process RNASEQ_GATK_RECALIBRATE { 
  tag "${sample_ID}"
  label "recab"

  input: 
    path genome
    path index
    path dict
    tuple val(sample_ID), path(bam), path(index)
    path(vcf_files)
    path(vcf_index_files)


  output:
    tuple val(sample_ID), path("${sample_ID}_final.uniq.bam"), path("${sample_ID}_final.uniq.bam.bai")

  script: 
  def variants = vcf_files.collect{ "$it" }.join(' --known-sites ')
  """
  gatk BaseRecalibrator \
          -R ${genome} \
          -I ${bam} \
          --known-sites ${variants} \
          -O final.rnaseq.grp 

  gatk ApplyBQSR \
          -R ${genome} -I ${bam} \
          --bqsr-recal-file final.rnaseq.grp \
          -O ${sample_ID}_final.uniq.bam

  samtools index ${sample_ID}_final.uniq.bam
  """
}


