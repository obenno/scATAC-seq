process SIGNAC {
    tag "${sample.id}"
    label 'process_high'
    publishDir "${params.outdir}/${sample.id}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(sample), path(fragment), path(fragmentIndex)
    path(gtf)

    output:
    tuple val(sample), path("${sample.id}_raw_cells.tsv.gz"),       emit: raw_cells
    tuple val(sample), path("${sample.id}_raw_meta.tsv.gz"),        emit: raw_meta
    tuple val(sample), path("${sample.id}_filtered_obj.rds"),       emit: obj
    tuple val(sample), path("${sample.id}_macs2_peaks.narrowPeak"), emit: macs_peaks

    script:
    def prefix = "${sample.id}"
    def fastqcThreads = Math.min(4, task.cpus)
    def nMem = "${task.memory.toBytes()}"
    def genomeSizeOpt = params.genomeSize ? "--genomeSize ${params.genomeSize}" : ""
    def blackListOpt = params.blackList ? "--blacklistBED ${params.blackList}" : ""
    def topCells = "${sample.expected_cells}"
    """
    signac_process.R -f $fragment \\
                     -g ${params.refGenome} \\
                     --gtf $gtf \\
                     -t $task.cpus \\
                     -m $nMem \\
                     ${genomeSizeOpt} \\
                     ${blackListOpt} \\
                     --emptyDrops_fdr ${params.emptyDrops_fdr} \\
                     --topCells ${topCells} \\
                     --minCell ${params.minCell} \\
                     --minFeature ${params.minFeature} \\
                     --nCount_min ${params.nCountMin} \\
                     --raw_cells_out ${sample.id}_raw_cells.tsv.gz \\
                     --raw_meta_metrics ${sample.id}_raw_meta.tsv.gz \\
                     --obj_out ${sample.id}_filtered_obj.rds
    mv macs2_peaks.narrowPeak ${sample.id}_macs2_peaks.narrowPeak
    """
}
