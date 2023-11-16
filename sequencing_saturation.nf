process CHECK_SATURATION {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(inputBam)
    tuple val(sampleID), path(inputBamBai)

    output:
    tuple val(sampleID), path("${sampleID}.saturation_out.json"), emit: outJSON

    script:
    // sinto fragments by default only recognize "chr" pattern in bam file
    if(params.refGenome == "hg38" || params.refGenome == "mm10"){
        chr_pattern = '(?i)^chr'
    }else if(params.refGenome == "hg38-mm10"){
        chr_pattern = '(?i)(^GRCh38_chr|^mm10___chr)'
    }else{
        chr_pattern = params.chrPattern
    }
    """
    get_sequencing_saturation.sh $inputBam "$chr_pattern" $task.cpus ${sampleID}.saturation_out.json
    """
}