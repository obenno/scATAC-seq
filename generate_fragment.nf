process GENERATE_FRAGMENTS {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/sinto/${sampleID}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    input:
    tuple val(sampleID), val(sampleName), path(bamFile)
    tuple val(sampleID), val(sampleName), path(bamIndex)
    //path whitelist

    output:
    tuple val(sampleID), val(sampleName), path("${sampleID}.fragments.sorted.bed.gz")       , emit: fragmentFile
    tuple val(sampleID), val(sampleName), path("${sampleID}.fragments.sorted.bed.gz.tbi")   , emit: fragmentIndex
    //path "versions.yml"                       , emit: versions

    script:
    """
    sinto fragments -b $bamFile -p $task.cpus -f ${sampleID}.fragments.bed --barcode_regex "[^:]*"
    sort -k1,1 -k2,2n ${sampleID}.fragments.bed > ${sampleID}.fragments.sorted.bed
    bgzip -@ $task.cpus ${sampleID}.fragments.sorted.bed
    tabix -p bed ${sampleID}.fragments.sorted.bed.gz
    rm ${sampleID}.fragments.bed
    """
}