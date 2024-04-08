process DEDUP {
    tag "${sampleID}"
    label 'process_low'

    input:
    tuple val(sampleID), path(inputBam), path(inputBamBai)

    output:
    tuple val(sampleID), path("${sampleID}.mark_dup.metrics"),  emit: dedup_metrics
    tuple val(sampleID), path("${sampleID}.mapping_stats.tsv"), emit: bam_stats

    script:
    """
    ## Mark duplicates
    export TMPDIR="./"
    ## Use samtools markdup instead of picard
    ## --no-multi-dup speed up the process
    samtools markdup -@ $task.cpus \\
                     -s \\
                     -f ${sampleID}.mark_dup.metrics \\
                     -m t \\
                     -t \\
                     --barcode-tag CB \\
                     --no-multi-dup \\
                     ${inputBam} \\
                     ${sampleID}.mark_dup.bam
    ## Generate bam stats
    samtools stats -@ $task.cpus ${sampleID}.mark_dup.bam > ${sampleID}.mapping_stats.tsv
 
    ## remove sorted bam to save space
    rm ${sampleID}.mark_dup.bam
    """
}
