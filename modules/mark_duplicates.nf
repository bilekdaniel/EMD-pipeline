process MARK_DUPLICATES {
  tag "${sample_ID}"
  label "mark_dup"

 input: 
    tuple val(sample_ID), path(bam), path(index)

  output:
    tuple val(sample_ID), path("${prefix}"), path(index)

  script: 

    prefix = task.ext.prefix ?: "${sample_ID}"
    def input_list = bam.collect{"--input $it"}.join(' ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}G -XX:-UsePerfData" \
        MarkDuplicatesSpark \
        $input_list \
        --output ${prefix} \
        --spark-master local[${task.cpus}] \
        --tmp-dir . \
    """
}

// at least 16 GB 
// max 16 cores - then diminishing returns
// 8 cores with 16 GB RAM might be good
// it is possible to collect metrics at the price of slow performance:  -M marked_dup_metrics.txt
// path("marked_dup_metrics.txt"), emit: metrics optional: true












// process MARK_DUPLICATES {
//   tag "${sample_ID}"
//   label "mark_dup"

//  input: 
//     tuple val(sample_ID), path(bam), path(index)

//   output:
//     tuple val(sample_ID), path("${sample_ID}_final.uniq.marked.bam"), path(index)

//   script: 
//   """
//   gatk MarkDuplicatesSpark \
//         -I ${bam} \
//         -O ${sample_ID}_final.uniq.marked.bam \
//         --tmp-dir . \
//         --spark-master local[${task.cpus}] 
//   """
// }

// // at least 16 GB 
// // max 16 cores - then diminishing returns
// // 8 cores with 16 GB RAM might be good
// // it is possible to collect metrics at the price of slow performance:  -M marked_dup_metrics.txt
// // path("marked_dup_metrics.txt"), emit: metrics optional: true
