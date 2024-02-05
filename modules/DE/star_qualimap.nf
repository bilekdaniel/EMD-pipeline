process STAR_QUALIMAP {
    label: "low"
    publishDir(path: "${params.output}${sample_ID}/qualimap", mode: 'copy', overwrite: 'true')

    input:
    tuple val(sample_ID), path(rgBAM)

    output:
    tuple val(sample_ID), path ("./qualimapReport.html"), emit: star_qualimap
    path "./css"
    path "./images_qualimapReport"
    path "./raw_data_qualimapReport"
    path "./rnaseq_qc_results.txt"

    """
    qualimap rnaseq -outdir ./ -bam ${rgBAM} -p strand-specific-reverse \
    -gtf ${params.gtf} --java-mem-size=10G
    """
}