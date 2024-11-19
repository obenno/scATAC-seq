process SIGNAC {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/${sampleID}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(sampleID), path(fragment), path(fragmentIndex)
    path(gtf)

    output:
    tuple val(sampleID), path("${sampleID}_raw_cells.tsv.gz"),    emit: raw_cells
    tuple val(sampleID), path("${sampleID}_raw_meta.tsv.gz"),    emit: raw_meta
    tuple val(sampleID), path("${sampleID}_filtered_obj.rds"), emit: obj
    tuple val(sampleID), path("${sampleID}_macs2_peaks.narrowPeak"), emit: macs_peaks

    script:
    def prefix = "${sampleID}"
    def fastqcThreads = Math.min(4, task.cpus)
    def nMem = "${task.memory.toBytes()}"
    def genomeSizeOpt = params.genomeSize ? "--genomeSize ${params.genomeSize}" : ""
    def blackListOpt = params.blackList ? "--blacklistBED ${params.blackList}" : ""
    """
    signac_process.R -f $fragment \\
                     -g ${params.refGenome} \\
                     --gtf $gtf \\
                     -t $task.cpus \\
                     -m $nMem \\
                     ${genomeSizeOpt} \\
                     ${blackListOpt} \\
                     --emptyDrops_fdr ${params.emptyDrops_fdr} \\
                     --minCell ${params.minCell} \\
                     --minFeature ${params.minFeature} \\
                     --nCount_min ${params.nCountMin} \\
                     --raw_cells_out ${sampleID}_raw_cells.tsv.gz \\
                     --raw_meta_metrics ${sampleID}_raw_meta.tsv.gz \\
                     --obj_out ${sampleID}_filtered_obj.rds
    mv macs2_peaks.narrowPeak ${sampleID}_macs2_peaks.narrowPeak
    """
}
