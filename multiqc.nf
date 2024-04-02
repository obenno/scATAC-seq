process MULTIQC {
    tag "${sampleID}"
    label 'process_low'
    publishDir "${params.outdir}/${sampleID}/multiqc",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/report.html/){
            return filename
        }else{
            return null
        }
    }

    input:
    tuple val(sampleID),
          path(read1_before),
          path(read2_before),
          path(read1_after),
          path(read2_after),
          path(cutadapt_report),
          path(bam_stats)

    output:
    tuple val(sampleID), path("${sampleID}_multiqc_report.html"), emit: report    

    script:
    def prefix = "${sampleID}"
    def fastqcThreads = Math.min(4, task.cpus)
    """
    ## rename input files
    ln -s ${read1_before} ${prefix}_raw_R1.fq.gz
    ln -s ${read2_before} ${prefix}_raw_R2.fq.gz
    ln -s ${read1_after} ${prefix}_trimmed_R1.fq.gz
    ln -s ${read2_after} ${prefix}_trimmed_R2.fq.gz
    fastqc -t $fastqcThreads --nogroup ${prefix}_raw_R1.fq.gz ${prefix}_raw_R2.fq.gz ${prefix}_trimmed_R1.fq.gz ${prefix}_trimmed_R2.fq.gz
    multiqc --title ${prefix} .
    """
}
