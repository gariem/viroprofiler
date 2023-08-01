
include { SPADES } from '../../modules/nf-core/spades/main'
include { MEGAHIT } from '../../modules/nf-core/megahit/main'


workflow ASSEMBLE {

    take:
        ch_clean_reads 

    main:

    // assembler options: "spades", "megahit"
    switch (params.assembler) {
        
        case "spades":
            SPADES (
                ch_clean_reads.map { meta, fastq -> [ meta, fastq, [], [] ] },
                []
            )
            contigs = SPADES.out.contigs
            scaffolds = SPADES.out.scaffolds
            versions = SPADES.out.versions.first()
            break;

        case "megahit":
            MEGAHIT ( 
                ch_clean_reads.map { meta, fastq -> [ meta, fastq ] }
            )
            contigs = MEGAHIT.out.contigs
            scaffolds = MEGAHIT.out.contigs
            versions = MEGAHIT.out.versions.first()

            break;
        
        default:
            exit 1, "Invalid assembler option: ${params.assembler}. Please specify one of [spades, megahit]"
            break;
    }

    // Use contigs or scaffolds
    if ( params.assemblies == "contigs") {
        assemblies = contigs
    } else {
        assemblies = scaffolds
    }

    emit:
        assemblies = assemblies
        versions = versions

}



