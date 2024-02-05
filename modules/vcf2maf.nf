
process VCF2MAF {
    tag "${sample_ID}"
    label "vcf2maf"
    publishDir(path: "${params.output}${sample_ID}", mode: 'copy', overwrite: 'true') 

    input: 
    tuple val(sample_ID), path(vcf), path(vcf_index)

    output:
    tuple val(sample_ID), path("${sample_ID}.vep.maf")

    script:
    """
    zcat ${vcf} > unziped.vcf
    perl /scratch/project/open-27-18/mrd_trans/genome/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf unziped.vcf --output-maf ${sample_ID}.vep.maf --vep-path /scratch/project/open-27-18/mrd_trans/work/conda/vcftomaf2-52140f1afe2fb579b114076c21d7e39b/bin --vep-data /scratch/project/open-27-18/mrd_trans/genome/.vep \
    --ref-fasta /scratch/project/open-27-18/mrd_trans/genome/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --species homo_sapiens --ncbi-build GRCh38 --cache-version 102 --tumor-id ${sample_ID}
    """
}

