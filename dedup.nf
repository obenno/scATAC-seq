process DEDUP {
    tag "${sample.id}"
    label 'process_low'

    input:
    tuple val(sample), path(inputBam), path(inputBamBai)

    output:
    tuple val(sample), path("${sample.id}.mark_dup.metrics"),  emit: dedup_metrics
    tuple val(sample), path("${sample.id}.mapping_stats.tsv"), emit: bam_stats

    script:
    """
    ## Mark duplicates
    export TMPDIR="./"
    ## Use samtools markdup instead of picard
    ## --no-multi-dup speed up the process
    samtools markdup -@ $task.cpus \\
                     -s \\
                     -f ${sample.id}.mark_dup.metrics \\
                     -m t \\
                     -t \\
                     --barcode-tag CB \\
                     --no-multi-dup \\
                     ${inputBam} \\
                     ${sample.id}.mark_dup.bam
    ## Generate bam stats
    samtools stats -@ $task.cpus ${sample.id}.mark_dup.bam > ${sample.id}.mapping_stats.tsv
 
    ## remove sorted bam to save space
    rm ${sample.id}.mark_dup.bam
    """
}
