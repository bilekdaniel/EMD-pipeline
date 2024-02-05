process PREPARE_VCF_FOR_ASE {
  label "low"
  tag "${sample_ID}"
  publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true') 

  input: 
    tuple val(sample_ID), path('final.vcf'), path('commonSNPs.diff.sites_in_files')
  output: 
    tuple val(sample_ID), path('known_snps.vcf.gz'), path('known_snps.vcf.gz.tbi')
    path 'AF.histogram.pdf'

  script:
  '''
  awk 'BEGIN{OFS="\t"} $4~/B/{print $1,$2,$3}' commonSNPs.diff.sites_in_files  > test.bed
    
  vcftools --vcf final.vcf --bed test.bed --recode --keep-INFO-all --stdout > known_snps.vcf

  grep -v '#'  known_snps.vcf | awk -F '\\t' '{print $10}' \
               |awk -F ':' '{print $2}'|perl -ne 'chomp($_); \
               @v=split(/\\,/,$_); if($v[0]!=0 ||$v[1] !=0)\
               {print  $v[1]/($v[1]+$v[0])."\\n"; }' |awk '$1!=1' \
               >AF.4R

  /scratch/project/open-27-18/mrd_trans/scripts/gghist.R -i AF.4R -o AF.histogram.pdf
  # Known SNPs have to be zipped and indexed for being used
  bgzip -c known_snps.vcf  > known_snps.vcf.gz
  tabix -p vcf known_snps.vcf.gz
  '''
}