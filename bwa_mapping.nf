process BWA_MAPPING {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read1), path(read2)
    path index

    output:
    tuple val(sampleID), path("${sampleID}.bwa_mapping.sorted.bam")           , emit: bam
    tuple val(sampleID), path("${sampleID}.bwa_mapping.sorted.bam.bai")       , emit: bai

    script:
    def prefix     = "${sampleID}"
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
    tag "${sampleID}"
    label 'process_low'

    input:
    tuple val(sampleID), path(inputBam), path(inputBamBai)

    output:
    tuple val(sampleID), path("${sampleID}.mapping_stats.tsv") ,emit: tsv

    script:
    def stats_threads = Math.min(2, task.cpus)
    """
    samtools stats -@ $task.cpus $inputBam > ${sampleID}.mapping_stats.tsv
    """
}