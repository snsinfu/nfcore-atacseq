process BEDTOOLS_GENOMECOV {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.0 bioconda::samtools=1.17 conda-forge::coreutils=9.3 conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0':
        'biocontainers/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' }"

    input:
    tuple val(meta), path(bam), path(flagstat)
    path(sizes)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    tuple val(meta), path("*.txt")     , emit: scale_factor
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def buffer = task.memory.toGiga().intdiv(2)
    """
    SCALE_FACTOR=\$(grep '[0-9] mapped (' $flagstat | awk '{print 1000000/\$1}')
    echo \$SCALE_FACTOR > ${prefix}.scale_factor.txt

    if [ "${meta.single_end}" = "true" ]; then
        bedtools \\
            genomecov \\
            -ibam $bam \\
            -bg \\
            -scale \$SCALE_FACTOR \\
            $args \\
        > tmp.bg
    else
        ## bedtools genomecov -pc silently drops pairs whose reverse mate starts
        ## upstream of the forward mate (dovetailing pairs, i.e. fragment length
        ## shorter than the read length after --ATACshift). Build one fragment
        ## [POS-1, POS-1+TLEN) per pair from its positive-TLEN (leftmost) mate.
        samtools view -f 0x2 -F 0x904 $bam \\
        | awk 'BEGIN{OFS="\\t"} \$9 > 0 {print \$3, \$4 - 1, \$4 - 1 + \$9}' \\
        | bedtools \\
            genomecov \\
            -i - \\
            -g $sizes \\
            -bg \\
            -scale \$SCALE_FACTOR \\
            $args \\
        > tmp.bg
    fi

    ## ref: https://www.biostars.org/p/66927/
    ## ref in nf-core: https://github.com/nf-core/hicar/blob/d2d17a924e42d6f88640b79d48d8b332f33a953f/modules/local/atacreads/bedsort.nf#L23-L29
    LC_ALL=C sort \\
        --parallel=$task.cpus \\
        --buffer-size=${buffer}G \\
        -k1,1 -k2,2n \\
        $args2 \\
        tmp.bg \\
    > ${prefix}.bedGraph

    rm tmp.bg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(samtools --version | head -n 1 | sed -e "s/samtools //g")
        sort: \$(sort --version | head -n 1 | awk '{print \$4;}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedGraph
    touch ${prefix}.scale_factor.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(samtools --version | head -n 1 | sed -e "s/samtools //g")
    END_VERSIONS
    """
}
