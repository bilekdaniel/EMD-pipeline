process STAR {
    label "star"
    publishDir(path: "${params.output}${sample_ID}/star", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), path(reads_wo_rrna1)
    tuple val(sample_ID), path(reads_wo_rrna2)

        
    output:
    tuple val(sample_ID), path("${sample_ID}_Aligned.sortedByCoord.out.bam"), path("${sample_ID}_Aligned.sortedByCoord.out.bam.bai") emit: BAM_FILE
    path "${sample_ID}_Log.final.out", optional: true
    path "${sample_ID}_Log.out", optional: true
    path "${sample_ID}_Log.progress.out", optional: true
    path "${sample_ID}_ReadsPerGene.out.tab", optional: true
    path "${sample_ID}_rgAdded_aligned.sortedByCoord.out", optional: true
    path "${sample_ID}_SJ.out.tab", optional: true

    """
    STAR --genomeDir ${params.star_index} \
    --runThreadN ${task.cpus} \
    --readFilesIn ${reads_wo_rrna1} ${reads_wo_rrna2} \
    --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
    --outSAMattributes Standard --quantMode GeneCounts \
    ${params.fasta_suffix_is_gz} --outFileNamePrefix ${sample_ID}_ \
    --chimSegmentMin 12 --chimJunctionOverhangMin 8 \
    --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 \
    --alignMatesGapMax 100000 --alignIntronMax 100000 \
    --chimOutType WithinBAM

    samtools index ${sample_ID}_Aligned.sortedByCoord.out.bam
    """
}

process RNASEQ_MAPPING_STAR {
  tag "${sample_ID}"
  label "star"

  input: 
    tuple val(sample_ID), path(reads_wo_rrna1)
    tuple val(sample_ID), path(reads_wo_rrna2) 
    path genome
    path genomeDir


  output: 
    tuple val(sample_ID), path('Aligned.sortedByCoord.uniq.bam'), path('Aligned.sortedByCoord.uniq.bam.bai')

  script:
  """
  # ngs-nf-dev Align reads to genome
  STAR --genomeDir ${genomeDir} \
       --readFilesIn ${reads_wo_rrna1} ${reads_wo_rrna2} \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999

  # Run 2-pass mapping (improve alignmets using table of splice junctions and create a new index)  
  STAR --genomeDir ${genomeDir} \
       --readFilesIn ${reads_wo_rrna1} ${reads_wo_rrna2} \
       --runThreadN ${task.cpus} \
       --readFilesCommand zcat \
       --outFilterType BySJout \
       --alignSJoverhangMin 8 \
       --alignSJDBoverhangMin 1 \
       --outFilterMismatchNmax 999 \
       --sjdbFileChrStartEnd SJ.out.tab \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMattrRGline ID:${sample_ID} LB:library PL:illumina PU:machine SM:GM12878

  # Select only unique alignments, no multimaps
  (samtools view -H Aligned.sortedByCoord.out.bam; samtools view Aligned.sortedByCoord.out.bam| grep -w 'NH:i:1') \
  |samtools view -Sb - > Aligned.sortedByCoord.uniq.bam
  
  # Index the BAM file
  samtools index Aligned.sortedByCoord.uniq.bam
  """
}