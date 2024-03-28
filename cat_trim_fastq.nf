process CHECK_BARCODE {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read1_list), path(read2_list)
    path(whitelist)

    output:
    tuple val(sampleID), path("${sampleID}_1.cutadapt_input.fq.gz"), emit: read1
    tuple val(sampleID), path("${sampleID}_2.cutadapt_input.fq.gz"), emit: read2
    tuple val(sampleID), path("${sampleID}_1.merged.fq.gz"),         emit: merged_r1
    tuple val(sampleID), path("${sampleID}_2.merged.fq.gz"),         emit: merged_r2
    tuple val(sampleID), path("${sampleID}_barcode_stats.tsv"),      emit: stats

    script:
    def prefix   = "${sampleID}"
    def read1 = read1_list.collect{ it.toString() }
    def read2 = read2_list.collect{ it.toString() }
    def whitelist_collapsed = whitelist.collect{ it.toString() }.join(',')
    def pigzThreads = Math.min(6, task.cpus)

    scriptString = []
    
    if(read1.size == 1 && read2.size == 1){
        scriptString.push(
        """
        ln -s ${read1.sort().join(' ')} ${prefix}_1.merged.fq.gz
        ln -s ${read2.sort().join(' ')} ${prefix}_2.merged.fq.gz
        """.stripIndent()
        )
    }else{
        scriptString.push(
        """
        cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fq.gz
        cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fq.gz
        """.stripIndent()
        )
    }

    scriptString.push(
    """
    extract_and_correct_thunderbio_barcode_from_fastq \
    $whitelist_collapsed \
    ${prefix}_1.merged.fq.gz \
    ${prefix}_2.merged.fq.gz \
    ${prefix}_1.barcode.fq \
    ${prefix}_2.barcode.fq \
    ${prefix}_barcode_stats.tsv \
    $task.cpus

    ## remove temp fastq
    rm combined_r1_r2_*.fq.gz
    
    ## compress fastq
    pigz -p $pigzThreads ${prefix}_1.barcode.fq ${prefix}_2.barcode.fq
    
    zcat ${prefix}_1.barcode.fq.gz | awk 'NR%4==1{if(\$4!="CB:Z:-"){print;getline;print;getline;print;getline;print}}' | pigz -p $pigzThreads > ${prefix}_1.cutadapt_input.fq.gz
    zcat ${prefix}_2.barcode.fq.gz | awk 'NR%4==1{if(\$4!="CB:Z:-"){print;getline;print;getline;print;getline;print}}' | pigz -p $pigzThreads > ${prefix}_2.cutadapt_input.fq.gz

    rm ${prefix}_1.barcode.fq.gz ${prefix}_2.barcode.fq.gz
    """.stripIndent()
    )

    scriptString.reverse().join()
}

process TRIM_FASTQ {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read1), path(read2)

    output:
    tuple val(sampleID),   path("*R1.trimmed.fq.gz"), emit: read1
    tuple val(sampleID),   path("*R2.trimmed.fq.gz"), emit: read2
    tuple val(sampleID),   path("${sampleID}.cutadapt.json"), emit: report_JSON

    script:
    def prefix   = "${sampleID}"
    """
    ## cutadapt QC and trim ME adapter
    cutadapt -j $task.cpus \\
             -q 30 \\
             -m $params.trim_mLen \\
             --json ${prefix}.cutadapt.json \\
             $params.trimOpt \\
             -o ${prefix}.R1.trimmed.fq.gz -p ${prefix}.R2.trimmed.fq.gz \\
             $read1 $read2
    """
}
