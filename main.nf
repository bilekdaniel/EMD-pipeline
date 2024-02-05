include { TRIM_GALORE_PE                                } from './modules/trim_galore'
include { RIBO_DET                                      } from './modules/ribo_det'
include { STAR                                          } from './modules/star'
include { MARK_DUPLICATES                               } from './modules/mark_duplicates.nf'
include { RNASEQ_GATK_SPLITNCIGAR                       } from './modules/gatk_splitncigar'
include { RNASEQ_GATK_RECALIBRATE                       } from './modules/gatk_recalibrate'
include { PREPARE_VCF_META                              } from './modules/prepare_vcf'
include { PREPARE_VCF_FILE; vcf_files_ch                } from './modules/prepare_vcf'
include { PREPARE_GENOME_SAMTOOLS                       } from './modules/prep_indexes'
include { PREPARE_GENOME_PICARD                         } from './modules/prep_indexes'
include { PREPARE_STAR_GENOME_INDEX                     } from './modules/prep_indexes'
include { RNASEQ_MAPPING_STAR                           } from './modules/star'
include { RNASEQ_CALL_VARIANTS                          } from './modules/rna_seq_variants'
include { POST_PROCESS_VCF                              } from './modules/post_process_vcf'
include { PREPARE_VCF_FOR_ASE                           } from './modules/prepare_vcf_for_ase'
include { ASE_KNOWNSNPS; group_per_sample               } from './modules/ase_knownsnp.nf'
include { VCF2MAF                                       } from './modules/vcf2maf'
include { MUTECT2                                       } from './modules/mutect2.nf'
include { MUTECT2_FILTER                                } from './modules/mutect2_filter.nf'
include { FUNCOTATOR                                    } from './modules/funcotator.nf'


Channel
    .fromPath(params.sampleCSV)
    .splitCsv(header: true)
    .map { row -> [sample_ID: row.names, fq1: "${params.samples}${row.fq1}", fq2: "${params.samples}${row.fq2}"] }
    .set { inv_samples }


workflow {

    PREPARE_GENOME_SAMTOOLS(params.genome)

    PREPARE_GENOME_PICARD(params.genome)

    PREPARE_STAR_GENOME_INDEX(params.genome)

    PREPARE_VCF_META(params.baserecalibrator_resources)

    PREPARE_VCF_FILE(vcf_files_ch (PREPARE_VCF_META.out) )

    TRIM_GALORE_PE(inv_samples)

    // RIBO_DET(TRIM_GALORE_PE.out)
    //  RNASEQ_MAPPING_STAR(RIBO_DET.out.reads_wo_rrna1, RIBO_DET.out.reads_wo_rrna2, params.genome, PREPARE_STAR_GENOME_INDEX.out)

    RNASEQ_MAPPING_STAR(TRIM_GALORE_PE.out, params.genome, PREPARE_STAR_GENOME_INDEX.out) //skipping ribodet

    MARK_DUPLICATES(RNASEQ_MAPPING_STAR.out)

    RNASEQ_GATK_SPLITNCIGAR(params.genome, PREPARE_GENOME_SAMTOOLS.out, PREPARE_GENOME_PICARD.out, MARK_DUPLICATES.out)

    RNASEQ_GATK_RECALIBRATE(params.genome, PREPARE_GENOME_SAMTOOLS.out, PREPARE_GENOME_PICARD.out, RNASEQ_GATK_SPLITNCIGAR.out, PREPARE_VCF_FILE.out.vcf.collect(), PREPARE_VCF_FILE.out.index.collect())

    MUTECT2(params.genome, PREPARE_GENOME_SAMTOOLS.out, PREPARE_GENOME_PICARD.out, params.POM, params.POM_INDEX, params.germline_resources, params.germline_resources_index, RNASEQ_GATK_RECALIBRATE.out)

    MUTECT2_FILTER(params.genome, PREPARE_GENOME_SAMTOOLS.out, PREPARE_GENOME_PICARD.out, MUTECT2.out)

    FUNCOTATOR(params.genome, PREPARE_GENOME_SAMTOOLS.out, PREPARE_GENOME_PICARD.out, params.funcotator_resource_folder, MUTECT2_FILTER.out )

//    VCF2MAF(MUTECT2_FILTER.out)
//    RNASEQ_CALL_VARIANTS(params.genome, PREPARE_GENOME_SAMTOOLS.out, PREPARE_GENOME_PICARD.out, RNASEQ_GATK_RECALIBRATE.out)
//
//    POST_PROCESS_VCF(RNASEQ_CALL_VARIANTS.out, PREPARE_VCF_FILE.out )

//    PREPARE_VCF_FOR_ASE(POST_PROCESS_VCF.out)

//    ASE_KNOWNSNPS(params.genome, PREPARE_GENOME_SAMTOOLS.out, PREPARE_GENOME_PICARD.out, group_per_sample(RNASEQ_GATK_RECALIBRATE.out, PREPARE_VCF_FOR_ASE.out[0]) )

//    VCF2MAF(RNASEQ_CALL_VARIANTS.out)


//    ADD_READ_GROUP(STAR.out.BAM_FILE)

//     STAR_QUALIMAP(ADD_READ_GROUP.out.rgBAM_FILE)

//     SALMON(RIBO_DET.out.reads_wo_rrna1, RIBO_DET.out.reads_wo_rrna2)


//    RNASEQ_GATK_RECALIBRATE(RNASEQ_GATK_SPLITNCIGAR.out.BAM, PREPARE_VCF_FILE.out.vcf)

//     multiqc (star_qualimap.out.star_qualimap.collect(), ribo_det.out.reads_wo_rrna1.collect(), ribo_det.out.reads_wo_rrna2.collect(), salmon.out.quant_folder.collect())

//   fastqc_both(ribo_det.out.reads_wo_rrna1, ribo_det.out.reads_wo_rrna2)

//    deseq2_paired(salmon.out.collect())

//    trinity (ribo_det.out.reads_wo_rrna1, ribo_det.out.reads_wo_rrna2)
}