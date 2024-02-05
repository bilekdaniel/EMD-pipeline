process ANNOTATE {
    tag "${sample_ID}"
    label "test"
    publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true') 
   
    input: 
    tuple val(sample_ID), path(vcf)

    output:
    tuple val(sample_ID), path("${sample_ID}.vep.vcf")

    script: 
    
"""
    vep --species homo_sapiens \
    --assembly GRCh38 \
    --offline \
    --no_progress \
    --no_stats \
    --sift b \
    --ccds \
    --uniprot \
    --hgvs \
    --symbol \
    --numbers \
    --domains \
    --gene_phenotype \
    --canonical \
    --protein \
    --biotype \
    --tsl \
    --pubmed \
    --variant_class \
    --shift_hgvs 1 \
    --check_existing \
    --total_length \
    --allele_number \
    --no_escape \
    --xref_refseq \
    --failed 1 \
    --vcf \
    --minimal \
    --flag_pick_allele \
    --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir ${params.genome}
    --fasta ${baseDir}/genome/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    --input_file ${vcf} \
    --output_file ${sample_ID}.vep.vcf \
    --polyphen b \
    --af \
    --af_1kg \
    --af_esp \
    --regulatory
    """
 }

//--dir ${baseDir}/genome/.vep \