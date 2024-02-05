process RNASEQ_CALL_VARIANTS {
  tag "${sample_ID}"
  label "medium"

  input:
    path genome
    path index
    path dict
    tuple val(sample_ID), path(bam), path(bai)
 
  output: 
    tuple val(sample_ID), path('final.vcf')

  script:

  """
  # fix absolute path in dict file
  sed -i 's@UR:file:.*${genome}@UR:file:${genome}@g' $dict
  
  # Variant calling
  gatk HaplotypeCaller \
          --native-pair-hmm-threads ${task.cpus} \
          --reference ${genome} \
          --output output.gatk.vcf.gz \
          -I ${bam} \
          --standard-min-confidence-threshold-for-calling 30.0 \
          --dont-use-soft-clipped-bases 

  # Variant filtering
  gatk VariantFiltration \
          -R ${genome} -V output.gatk.vcf.gz \
          --cluster-window-size 35 --cluster-size 3 \
          --filter-name FS --filter-expression \"FS > 30.0\" \
          --filter-name QD --filter-expression \"QD < 2.0\" \
          -O final.vcf
  """
}
