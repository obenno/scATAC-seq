process CAT_FASTQ {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read1_list), path(read2_list)

    output:
    tuple val(sampleID), path("${sampleID}_1.merged.fq.gz"),         emit: read1
    tuple val(sampleID), path("${sampleID}_2.merged.fq.gz"),         emit: read2

    script:
    def prefix   = "${sampleID}"
    def read1 = read1_list.collect{ it.toString() }
    def read2 = read2_list.collect{ it.toString() }
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
    scriptString.reverse().join()
    
}


process CHECK_BARCODE {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read1), path(read2)
    path(whitelist)

    output:
    tuple val(sampleID), path("${sampleID}_1.cutadapt_input.fq.gz"), emit: read1
    tuple val(sampleID), path("${sampleID}_2.cutadapt_input.fq.gz"), emit: read2
    tuple val(sampleID), path("${sampleID}_barcode_stats.tsv"),      emit: stats

    script:
    def prefix   = "${sampleID}"
    def whitelist_collapsed = whitelist.collect{ it.toString() }.join(',')
    def pigzThreads = Math.min(6, task.cpus)

    """
    $baseDir/codon/bin/codon run -release -plugin seq $baseDir/bin/extract_and_correct_thunderbio_barcode_from_fastq.codon \\
    $whitelist_collapsed \\
    ${read1} \\
    ${read2} \\
    ${prefix}_1.barcode.fq \\
    ${prefix}_2.barcode.fq \\
    ${prefix}_barcode_stats.tsv \\
    $task.cpus

    ## remove temp fastq
    rm combined_r1_r2_*.fq.gz
    
    awk 'NR%4==1{if(\$4!="CB:Z:-"){print;getline;print;getline;print;getline;print}}' ${prefix}_1.barcode.fq |
        pigz -p $pigzThreads > ${prefix}_1.cutadapt_input.fq.gz
    awk 'NR%4==1{if(\$4!="CB:Z:-"){print;getline;print;getline;print;getline;print}}' ${prefix}_2.barcode.fq |
        pigz -p $pigzThreads > ${prefix}_2.cutadapt_input.fq.gz

    rm ${prefix}_1.barcode.fq ${prefix}_2.barcode.fq
    """
}

process CAT_FASTQ_10X {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read1_list), path(read2_list), path(read3_list)

    output:
    tuple val(sampleID), path("${sampleID}_1.merged.fq.gz"),         emit: read1
    tuple val(sampleID), path("${sampleID}_2.merged.fq.gz"),         emit: read2
    tuple val(sampleID), path("${sampleID}_3.merged.fq.gz"),         emit: read3

    script:
    def prefix   = "${sampleID}"
    def read1 = read1_list.collect{ it.toString() }
    def read2 = read2_list.collect{ it.toString() }
    def read3 = read3_list.collect{ it.toString() }
    scriptString = []
    
    if(read1.size == 1 && read2.size == 1 && read3.size == 1){
        scriptString.push(
        """
        ln -s ${read1.sort().join(' ')} ${prefix}_1.merged.fq.gz
        ln -s ${read2.sort().join(' ')} ${prefix}_2.merged.fq.gz
        ln -s ${read3.sort().join(' ')} ${prefix}_3.merged.fq.gz
        """.stripIndent()
        )
    }else{
        scriptString.push(
        """
        cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fq.gz
        cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fq.gz
        cat ${read3.sort().join(' ')} > ${prefix}_3.merged.fq.gz
        """.stripIndent()
        )
    }
    scriptString.reverse().join()
    
}

process CHECK_BARCODE_10X {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read2), path(read1), path(read3)
    path(whitelist)

    output:
    tuple val(sampleID), path("${sampleID}_1.cutadapt_input.fq.gz"), emit: read1
    tuple val(sampleID), path("${sampleID}_2.cutadapt_input.fq.gz"), emit: read2
    tuple val(sampleID), path("${sampleID}_barcode_stats.tsv"),      emit: stats

    script:
    def prefix   = "${sampleID}"
    def pigzThreads = Math.min(6, task.cpus)

    """
    $baseDir/codon/bin/codon run -release -plugin seq $baseDir/bin/extract_and_correct_10x_barcode_from_fastq.codon \\
    ${whitelist} \\
    ${read2} \\
    ${read1} \\
    ${read3} \\
    ${prefix}_1.barcode.fq \\
    ${prefix}_2.barcode.fq \\
    ${prefix}_barcode_stats.tsv \\
    $task.cpus

    ## remove temp fastq
    rm combined_r1_r2_*.fq.gz
    
    awk 'NR%4==1{if(\$4!="CB:Z:-"){print;getline;print;getline;print;getline;print}}' ${prefix}_1.barcode.fq |
        pigz -p $pigzThreads > ${prefix}_1.cutadapt_input.fq.gz
    awk 'NR%4==1{if(\$4!="CB:Z:-"){print;getline;print;getline;print;getline;print}}' ${prefix}_2.barcode.fq |
        pigz -p $pigzThreads > ${prefix}_2.cutadapt_input.fq.gz

    rm ${prefix}_1.barcode.fq ${prefix}_2.barcode.fq
    """
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
