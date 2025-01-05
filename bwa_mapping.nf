process BWA_MAPPING {
    tag "${sample.id}"
    label 'process_high'

    input:
    tuple val(sample), path(read1), path(read2)
    path index

    output:
    tuple val(sample), path("${sample.id}.bwa_mapping.sorted.bam")           , emit: bam
    tuple val(sample), path("${sample.id}.bwa_mapping.sorted.bam.bai")       , emit: bai

    script:
    def prefix     = "${sample.id}"
    def genomeIndex = index.collect{ it.toString() }[0].split("\\.")[0]
    """
    bwa-mem2 mem \\
    -t $task.cpus \\
    -C \\
    $genomeIndex \\
    $read1 \\
    $read2 |
        samtools fixmate -u -m -O bam - - |
        samtools view -b -@ $task.cpus > ${prefix}.bwa_mapping.bam
    samtools sort -@ $task.cpus ${prefix}.bwa_mapping.bam > ${prefix}.bwa_mapping.sorted.bam
    rm ${prefix}.bwa_mapping.bam
    samtools index ${prefix}.bwa_mapping.sorted.bam
    """
}

process BAM_SUMMARY {
    tag "${sample.id}"
    label 'process_low'

    input:
    tuple val(sample), path(inputBam), path(inputBamBai)

    output:
    tuple val(sample), path("${sample.id}.mapping_stats.tsv") ,emit: tsv

    script:
    def stats_threads = Math.min(2, task.cpus)
    """
    samtools stats -@ $task.cpus $inputBam > ${sample.id}.mapping_stats.tsv
    """
}
