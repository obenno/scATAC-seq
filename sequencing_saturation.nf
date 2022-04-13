process CHECK_SATURATION {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), val(sampleName), path(inputBam)
    tuple val(sampleID), val(sampleName), path(inputBamBai)

    output:
    tuple val(sampleID), val(sampleName), path("${sampleID}.saturation_out.json"), emit: outJSON

    shell:
    """
    get_sequencing_saturation.sh $inputBAM $task.cpus ${sampleID}.saturation_out.json
    """
}