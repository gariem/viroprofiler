
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { FASTP } from '../../modules/nf-core/fastp/main'
include { DECONTAM } from '../../modules/local/decontam'

workflow PREPROCESS {

    take:
        reads 

    main:
        ch_versions = Channel.empty()
        ch_multiqc_files = Channel.empty()

        // MODULE: Run FastQC
        FASTQC (
            reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())

        // MODULE: Run fastp
        FASTP (
            reads, false, false 
        )
        ch_versions = ch_versions.mix(FASTP.out.versions.first())

        // Decontamination
        if (params.use_decontam) {
            ch_contamref = Channel.fromPath("${params.contamref_idx}", checkIfExists: true).first()
            DECONTAM (FASTP.out.reads, ch_contamref)
            ch_clean_reads = DECONTAM.out.reads
            ch_versions = ch_versions.mix(DECONTAM.out.versions.first())
        } else {
            ch_clean_reads = FASTP.out.reads
        }

        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))

    emit:
        ch_clean_reads = FASTP.out.reads
        versions = ch_versions
        ch_multiqc_files = ch_multiqc_files

}