process PREPARE_GENOME_SAMTOOLS { 
  tag "$genome.baseName"
  label "low"
  
  publishDir(path: "${params.genome_folder}", mode: 'copy', overwrite: 'true')
  input: 
    path genome
 
  output: 
    path "${genome}.fai"
  
  script:
  """
  samtools faidx ${genome}
  """
}


process PREPARE_GENOME_PICARD {
  tag "$genome.baseName"
  label "low"

  publishDir(path: "${params.genome_folder}", mode: 'copy', overwrite: 'true')
  input:
    path genome
  output:
    path "${genome.baseName}.dict"

  script:
  """
  gatk CreateSequenceDictionary -R ${genome} -O ${genome.baseName}.dict
  """
}


process PREPARE_STAR_GENOME_INDEX {
  tag "$genome.baseName"
  label "index_star"

  publishDir(path: "${params.genome_folder}", mode: 'copy', overwrite: 'true')

  input:
    path genome
  output:
    path "genome_dir"

  script:
  """
  mkdir -p genome_dir

  STAR --runMode genomeGenerate \
       --genomeDir genome_dir \
       --genomeFastaFiles ${genome} \
       --runThreadN ${task.cpus}
  """
}
//add sjdbGTFfiles