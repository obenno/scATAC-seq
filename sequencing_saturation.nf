process CHECK_SATURATION {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(cellsTSV), path(peaks), path(fragments), path(fragIndex)

    output:
    tuple val(sampleID), path("${sampleID}.saturation_out.json"), emit: outJSON

    script:
    """
    get_sequencing_saturation.v2.sh --frags ${fragments} --cells ${cellsTSV} --peaks ${peaks} --out ${sampleID}.saturation_out.json --threads $task.cpus
    """
}
