process DRIMSeq {
    label "high"
    publishDir(path: "${params.baseDir}/DRIMSeq", mode: 'copy', overwrite: 'true')

    input:
    val sample_ID
    path dir
    path csv
    path gtf
    path run_name
    val alpha
    val mge
    val mfe
    val mfp

    output:
    path("drimseq_completed")

    script:
    """
    Rscript bin/drimseq.R \
    ${dir} \ 
    ${csv} \
    ${gtf} \
    ${run_name} \ 
    ${alpha} \ 
    ${mge} \
    ${mfe} \
    ${mfp}

    echo drimseq_completed > drimseq_completed
    """
}
//Rscript bin/drimseq.R ./ data/samples_all.csv ../gencode/gencode.v36.annotation.gtf "" 0.05 10 10 0.1
workflow{
    DRIMSeq (
        salmon.out.map { it[0] }.collect(),
        params.baseDir,
        params.sample.sampleCSV,
        params.gtf,
        params.drimseq_run ?: "",
        params.drimseq_alpha ?: "0.05",
        params.drimseq_mge ?: "10",
        params.drimseqmfe ?: "10",
        params.drimseq_mfp ?: "0.1"
    )

}