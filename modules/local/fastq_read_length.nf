/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sample the first records of a FASTQ to estimate the read length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTQ_READ_LENGTH {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.read_length.txt"), emit: read_length

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "$fastq" == *.gz ]]; then CAT="zcat"; else CAT="cat"; fi
    # 'head' closes the pipe early; run without pipefail so the SIGPIPE (141) is not an error
    set +o pipefail
    \$CAT "$fastq" | head -n 40000 | awk 'NR % 4 == 2 { print length(\$0) }' | sort -nr | head -n 1 > ${prefix}.read_length.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo 150 > ${prefix}.read_length.txt
    """
}
