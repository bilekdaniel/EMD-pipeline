process RNASEQ_GATK_SPLITNCIGAR {
  tag "${sample_ID}"
  label 'medium'

  input: 
    path genome
    path index
    path genome_dict
    tuple val(sample_ID), path(bam), path(index)

  output:
    tuple val(sample_ID), path("${sample_ID}_split.bam"), path("${sample_ID}_split.bai")
  
  script:
  """
  gatk SplitNCigarReads \
            -R ${genome} \
            -I ${bam} \
            -O ${sample_ID}_split.bam
  """
}
// removed this argument: --refactor-cigar-string"