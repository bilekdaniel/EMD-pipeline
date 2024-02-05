process POST_PROCESS_VCF {
  tag "${sample_ID}"
  label "low"

  publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true')

  input:
    tuple val(sample_ID), path('final.vcf')
    tuple path('filtered.recode.vcf.gz'), path('filtered.recode.vcf.gz.tbi')
    
  output: 
    tuple val(sample_ID), path('final.vcf'), path('commonSNPs.diff.sites_in_files')
  
  script:
  '''

  grep -v '#' final.vcf | awk '$7~/PASS/' |perl -ne 'chomp($_); ($dp)=$_=~/DP\\=(\\d+)\\;/; if($dp>=8){print $_."\\n"};' > result.DP8.vcf

  awk '{print $1}' result.DP8.vcf | sort -u > unique.elements.txt

  vcftools --vcf result.DP8.vcf \
  --gzdiff filtered.recode.vcf.gz \
  --chr 1 \
  --chr 10 \
  --chr 11 \
  --chr 12 \
  --chr 13 \
  --chr 14 \
  --chr 15 \
  --chr 16 \
  --chr 17 \
  --chr 18 \
  --chr 19 \
  --chr 2 \
  --chr 20 \
  --chr 21 \
  --chr 22 \
  --chr 3 \
  --chr 4 \
  --chr 5 \
  --chr 6 \
  --chr 7 \
  --chr 8 \
  --chr 9 \
  --chr X \
  --diff-site --out commonSNPs
  '''
}