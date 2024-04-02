process SIGNAC {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/${sampleID}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/cells_out\.tsv|\.rds/){
            return filename
        }else{
            return null
        }
    }

    input:
    tuple val(sampleID), path(fragment), path(fragmentIndex)
    path(gtf)

    output:
    tuple val(sampleID), path("${sampleID}_raw_cells.tsv"),    emit: raw_cells
    tuple val(sampleID), path("${sampleID}_raw_meta.tsv"),    emit: raw_meta
    tuple val(sampleID), path("${sampleID}_filtered_obj.rds"), emit: obj
    tuple val(sampleID), path("${sampleID}_macs2_peaks.narrowPeak"), emit: macs_peaks

    script:
    def prefix = "${sampleID}"
    def fastqcThreads = Math.min(4, task.cpus)
    def nMem = "${task.memory.toBytes()}"
    """
    signac_process.R -f $fragment \\
                     -g ${params.refGenome} \\
                     --gtf $gtf \\
                     -t $task.cpus \\
                     -m $nMem \\
                     --nCount_min ${params.nCountMin} \\
                     --raw_cells_out ${sampleID}_raw_cells.tsv \\
                     --raw_meta_metrics ${sampleID}_raw_meta.tsv \\
                     --obj_out ${sampleID}_filtered_obj.rds
    mv macs2_peaks.narrowPeak ${sampleID}_macs2_peaks.narrowPeak
    """
}