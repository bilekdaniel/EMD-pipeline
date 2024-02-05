process ASE_KNOWNSNPS {
  tag "${sample_ID}"
  publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true') 
  label "medium"

  input:
    path genome
    path index
    path dict
    tuple val(sample_ID), path(vcf), path(tbi), path(bam), path(bai)
  
  output:
    path "ASE.tsv"
  
  script:
  """

  gatk SelectVariants -R ${genome} --variant ${vcf} --restrict-alleles-to BIALLELIC -select 'vc.getHetCount()==1' --select-type-to-include SNP -O sorted.vcf.gz

  gatk ASEReadCounter \
          -R ${genome} \
          -O ASE.tsv \
          -I ${bam} \
          -V sorted.vcf.gz
  """
}

def group_per_sample(bam_ch, vcf_ch) {
  bam_ch
    .join(vcf_ch)
    .view()
    .map{ joined -> 
      def sample_ID = joined[0]
      def bam = joined[1]
      def bai = joined[2]
      def vcf = joined[3]
      def tbi = joined[4]
      tuple(sample_ID, vcf, tbi, bam, bai)
    }
}
