process CAT_TRIM_FASTQ {
    tag "${sampleID}"
    label 'process_medium'
    publishDir "${params.outdir}/cutqc",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    // conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    // } else {
    //     container "biocontainers/biocontainers:v1.2.0_cv1"
    // }

    input:
    tuple val(sampleID), val(sampleName), path(read1_list), path(read2_list)
    path(whitelist)

    output:
    tuple val(sampleID), val(sampleName), path("*_1.merged.cDNA_R1.barcoded.fq.gz"), emit: read1
    tuple val(sampleID), val(sampleName), path("*_2.merged.barcoded.fq.gz"), emit: read2
    tuple val(sampleID), val(sampleName), path("*_1.merged.bc.fq.gz"), emit: read_bc
    tuple val(sampleID), val(sampleName), path("${sampleID}.readReport.json"), emit: readReport
    //tuple val(sampleID), val(sampleName), path("*.read_stat.txt"), emit: report
    //tuple val(sampleID), val(sampleName), path("*_cutqc_report.html"), emit: cutqc_report
    //path "versions.yml"                       , emit: versions

    script:
    def prefix   = "${sampleID}"
    def read1 = read1_list.collect{ it.toString() }
    def read2 = read2_list.collect{ it.toString() }
    """
    cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fq.gz
    cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fq.gz

    split_barcode_reads.sh ${prefix}_1.merged.fq.gz $params.bcReadLen

    ## use bash script instead, will be much faster but require more memory
    ##sinto barcode -b $params.bcLen --barcode_fastq ${prefix}_1.merged.bc.fq.gz \\
    ##--read1 ${prefix}_1.merged.cDNA_R1.fq.gz \\
    ##--read2 ${prefix}_2.merged.fq.gz
    add_barcode_to_reads.sh ${sampleID} ${whitelist} $params.bcLen ${prefix}_1.merged.bc.fq.gz ${prefix}_1.merged.cDNA_R1.fq.gz ${prefix}_2.merged.fq.gz > ${sampleID}.readReport.json

    ## rename sinto output to fq.gz
    ##mv ${prefix}_1.merged.cDNA_R1.barcoded.fastq.gz ${prefix}_1.merged.cDNA_R1.barcoded.fq.gz
    ##mv ${prefix}_2.merged.barcoded.fastq.gz ${prefix}_2.merged.barcoded.fq.gz

    ##cutqc cutqc ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz \\
    ##${prefix}_cutqc_report.html \\
    ##-j $task.cpus \\
    ##-m $params.trimLength \\
    ##-q $params.trimQuality \\
    ##$params.cutadaptOption

    ## remove merged fq.gz to save space
    rm ${prefix}_1.merged.fq.gz ${prefix}_2.merged.fq.gz
    """
}
