process GENERATE_FRAGMENTS {
    tag "${sample.id}"
    label 'process_high'
    publishDir "${params.outdir}/${sample.id}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
        
    input:
    tuple val(sample), path(bamFile), path(bamIndex)

    output:
    tuple val(sample), path("${sample.id}.fragments.sorted.bed.gz")       , emit: fragmentFile
    tuple val(sample), path("${sample.id}.fragments.sorted.bed.gz.tbi")   , emit: fragmentIndex

    script:
    // sinto fragments by default only recognize "chr" pattern in bam file
    if(params.refGenome == "hg38" || params.refGenome == "mm10"){
        chr_pattern = '(?i)^chr'
    }else if(params.refGenome == "hg38-mm10"){
        chr_pattern = '(?i)(^GRCh38_chr|^mm10___chr)'
    }else{
        chr_pattern = params.chrPattern
    }

    // use barcodetag or barcode_regex option to extract barcode
    def barcode_opts = params.barcode_regex ? "--barcode_regex $params.barcode_regex" : "--barcodetag $params.barcodetag"
    """
    export TMPDIR="./"
    ## Added --collapse_within option, see https://github.com/timoast/sinto/issues/36
    sinto fragments \\
    -b $bamFile \\
    -p $task.cpus \\
    -f ${sample.id}.fragments.bed \\
    ${barcode_opts} \\
    --use_chrom "${chr_pattern}" \\
    --collapse_within \\
    --max_distance ${params.max_distance} \\
    --min_distance ${params.min_distance} \\
    --chunksize ${params.chunksize}
    sort -k1,1 -k2,2n --parallel $task.cpus -S 256M ${sample.id}.fragments.bed > ${sample.id}.fragments.sorted.bed
    rm ${sample.id}.fragments.bed
    bgzip -@ $task.cpus ${sample.id}.fragments.sorted.bed
    tabix -p bed ${sample.id}.fragments.sorted.bed.gz
    """
}
