process BOWTIE2 {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/bowtie2/${sampleID}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    input:
    tuple val(sampleID), val(sampleName), path(read1)
    tuple val(sampleID), val(sampleName), path(read2)
    path index
    //path gtf
    //path whitelist

    output:
    tuple val(sampleID), val(sampleName), path("${sampleID}.sorted.bam")       , emit: bam
    tuple val(sampleID), val(sampleName), path("${sampleID}.sorted.bam.bai")   , emit: bai
    tuple val(sampleID), val(sampleName), path("${sampleID}.mappingReport.json")   , emit: mappingReport
    //path "versions.yml"                       , emit: versions

    script:
    def prefix     = "${sampleID}"
    def genomeIndex = index.collect{ it.toString() }[0].replaceFirst(/(\.rev)*\.(\d+)\.bt2(l)*$/, "")
    //def barcodeMate = params.bc_read == "fastq_1" ? 1 : 2
    // support pair-end only for now
    //sample.single_end = false
    """
    bowtie2 --mm \\
    -p $task.cpus \\
    -x $genomeIndex \\
    -X 2000 \\
    -1 $read1 \\
    -2 $read2 2> ${sampleID}.mappingReport.txt |
        samtools view -b |
        samtools sort -@ $task.cpus -o ${sampleID}.sorted.bam

    samtools index *.bam

    ## convert report.txt to json
    bowtie2_report_to_json.sh ${sampleID} ${sampleID}.mappingReport.txt > ${sampleID}.mappingReport.json
    """
}
