process DEDUP {
    tag "${sampleID}"
    label 'process_low'
    //publishDir "${params.outdir}/processedBam/${sampleID}",
    //    mode: "${params.publish_dir_mode}",
    //    enabled: params.outdir as boolean

    input:
    tuple val(sampleID), val(sampleName), path(inputBam)
    tuple val(sampleID), val(sampleName), path(inputBamBai)

    output:
    tuple val(sampleID), val(sampleName), env(ratio), emit: dupRatio
    //tuple val(sampleID), val(sampleName), path("${sampleID}.final.bam"), emit: finalBam
    //tuple val(sampleID), val(sampleName), path("${sampleID}.final.bam.bai"), emit: finalBai

    script:
    """
    ## Mark duplicates
    export TMPDIR="./"
    sortedBAM=\$(mktemp -p ./ sorted.XXXXXX.bam)
    samtools sort -@ $task.cpus ${inputBam} > \$sortedBAM
    mkdir picard_tmp
    picard MarkDuplicates \\
        TMP_DIR=\$(pwd)/picard_tmp \\
        INPUT=\$sortedBAM OUTPUT=${sampleID}.mark_dup.bam METRICS_FILE=${sampleID}.mark_dup.metrics \\
        VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=coordinate
    ## Filter alignments (remove duplicates)
    ##samtools view -@ $task.cpus -b -f 0x2 -F 0x4 -F 0x8 -F 0x100 -F 0x800 -F 0x400 -q 30 ${sampleID}.mark_dup.bam |
    ##    samtools sort -@ $task.cpus > ${sampleID}.final.bam
    ##samtools index ${sampleID}.final.bam
    ## remove sorted bam to save space
    rm \$sortedBAM
    rm ${sampleID}.mark_dup.bam
    ratio=\$(awk '\$0~/^## METRICS CLASS/{getline; getline; pos=NF-1; print \$pos}' ${sampleID}.mark_dup.metrics)
    """
}
