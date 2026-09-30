/*
 * Map reads with BWA, bwa-mem2 or bwa-mem3; sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { BWA_MEM                 } from '../../modules/nf-core/bwa/mem/main'
include { BWAMEM2_MEM             } from '../../modules/local/bwamem2_mem'
include { BWAMEM3_MEM             } from '../../modules/local/bwamem3_mem'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'

workflow ALIGN_BWA {
    take:
    ch_reads     // channel: [ val(meta), [ reads ] ]
    ch_index     // channel: /path/to/bwa|bwamem2|bwamem3/index/
    ch_fasta     // channel: /path/to/fasta
    val_aligner  //  string: 'bwa', 'bwa-mem2' or 'bwa-mem3'
    val_sort_bam // boolean: whether the aligner module sorts the output

    main:
    ch_versions = Channel.empty()
    ch_bam_orig = Channel.empty()

    if (val_aligner == 'bwa') {
        BWA_MEM ( ch_reads, ch_index, ch_fasta, val_sort_bam )
        ch_bam_orig = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    } else if (val_aligner == 'bwa-mem2') {
        BWAMEM2_MEM ( ch_reads, ch_index, ch_fasta, val_sort_bam )
        ch_bam_orig = BWAMEM2_MEM.out.bam
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    } else if (val_aligner == 'bwa-mem3') {
        BWAMEM3_MEM ( ch_reads, ch_index, ch_fasta, val_sort_bam )
        ch_bam_orig = BWAMEM3_MEM.out.bam
        ch_versions = ch_versions.mix(BWAMEM3_MEM.out.versions.first())
    } else {
        error "ALIGN_BWA: unsupported aligner '${val_aligner}'. Expected one of: bwa, bwa-mem2, bwa-mem3."
    }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( ch_bam_orig, ch_fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig = ch_bam_orig                          // channel: [ val(meta), bam      ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}
