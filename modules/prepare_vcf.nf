process PREPARE_VCF_META {

  input:
  path('resources_path')

  output:
  path("output.csv")

  script:
  '''
  folder_path=resources_path
  echo "ID,files" > output.csv

  index=1
  for file in "$folder_path"/*.gz; do
    absolute_path="$(realpath "$file")"
    echo "$index,$absolute_path" >> output.csv
    ((index++))
  done
  '''
}

def vcf_files_ch(resources_path) {
  resources_path
    .splitCsv(header: true)
    .map { row -> [file_ID: row.ID, files: row.files] }
}

process PREPARE_VCF_FILE {
label "low"
tag "$file_ID"
publishDir(path: "${params.baserecalibrator_resources}", mode: 'copy', overwrite: 'true')

  input: 
    tuple val(file_ID), path(files)
  
  output:
    path(files), emit: vcf
    path("${files}.tbi"), emit: index
  
  script:  
  """
  tabix ${files}
  """
}

