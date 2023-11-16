process GENERATE_FRAGMENTS {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/sinto/${sampleID}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    input:
    tuple val(sampleID), path(bamFile)
    tuple val(sampleID), path(bamIndex)
    //path whitelist

    output:
    tuple val(sampleID), path("${sampleID}.fragments.sorted.bed.gz")       , emit: fragmentFile
    tuple val(sampleID), path("${sampleID}.fragments.sorted.bed.gz.tbi")   , emit: fragmentIndex
    //path "versions.yml"                       , emit: versions

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
    export TMPDIR="./"
    sinto fragments -b $bamFile -p $task.cpus -f ${sampleID}.fragments.bed --barcode_regex "[^:]*" --use_chrom "${chr_pattern}"
    sort -k1,1 -k2,2n ${sampleID}.fragments.bed > ${sampleID}.fragments.sorted.bed
    bgzip -@ $task.cpus ${sampleID}.fragments.sorted.bed
    tabix -p bed ${sampleID}.fragments.sorted.bed.gz
    rm ${sampleID}.fragments.bed
    """
}
