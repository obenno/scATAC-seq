process CHECK_SATURATION {
    tag "${sample.id}"
    label 'process_high'
    publishDir "${params.outdir}/${sample.id}/saturation",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean

    input:
    tuple val(sample), path(cellsTSV), path(peaks), path(fragments), path(fragIndex)

    output:
    tuple val(sample), path("${sample.id}.saturation_out.json"), emit: outJSON

    script:
    """
    get_sequencing_saturation.v2.sh --frags ${fragments} --cells ${cellsTSV} --peaks ${peaks} --out ${sample.id}.saturation_out.json --threads $task.cpus
    """
}
