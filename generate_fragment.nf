process GENERATE_FRAGMENTS {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/${sampleID}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
        
    input:
    tuple val(sampleID), path(bamFile), path(bamIndex)

    output:
    tuple val(sampleID), path("${sampleID}.fragments.sorted.bed.gz")       , emit: fragmentFile
    tuple val(sampleID), path("${sampleID}.fragments.sorted.bed.gz.tbi")   , emit: fragmentIndex

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
    
    """
    export TMPDIR="./"
    ## Added --collapse_within option, see https://github.com/timoast/sinto/issues/36
    sinto fragments \\
    -b $bamFile \\
    -p $task.cpus \\
    -f ${sampleID}.fragments.bed \\
    --barcodetag ${params.barcodetag}  \\
    --barcode_regex ${params.barcode_regex} \\
    --use_chrom "${chr_pattern}" \\
    --collapse_within \\
    --max_distance ${params.max_distance} \\
    --min_distance ${params.min_distance} \\
    --chunksize ${params.chunksize}
    sort -k1,1 -k2,2n --parallel $task.cpus -S 256M ${sampleID}.fragments.bed > ${sampleID}.fragments.sorted.bed
    rm ${sampleID}.fragments.bed
    bgzip -@ $task.cpus ${sampleID}.fragments.sorted.bed
    tabix -p bed ${sampleID}.fragments.sorted.bed.gz
    """
}
