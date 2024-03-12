process CHECK_BARCODE {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(read1_list), path(read2_list)
    path(gtf)
    path(starIndex)
    path(whitelist)

    output:
    tuple val(sampleID),   path("${sampleID}_1.merged.fq.gz"), emit: merged_read1
    tuple val(sampleID),   path("${sampleID}_2.merged.fq.gz"), emit: merged_read2
    tuple val(sampleID),   path("*Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sampleID),   path("${sampleID}.Solo.out/Gene/Summary.csv"), emit: summary

    script:
    def prefix   = "${sampleID}"
    def read1 = read1_list.collect{ it.toString() }
    def read2 = read2_list.collect{ it.toString() }
    def pigzThreads = Math.min(6, task.cpus)
    
    scriptString = []

    scriptString.push(
    """
    ## set limit for concurrent opening files
    ulimit -Sn 40960
    """.stripIndent()
    )

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
    ## Perform barcode detection and correction by STARsolo
    STAR --runThreadN ${task.cpus} \\
    --genomeDir $starIndex \\
    --sjdbGTFfile $gtf \\
    --soloBarcodeMate 0 \\
    --readFilesIn ${prefix}_2.merged.fq.gz ${prefix}_1.merged.fq.gz \\
    --soloCBwhitelist $whitelist \\
    --soloBarcodeReadLength 0 \\
    --readFilesCommand zcat \\
    --outFileNamePrefix ${prefix}. \\
    --soloStrand $params.soloStrand \\
    --outFilterScoreMin 0 \\
    --soloCBmatchWLtype 1MM \\
    --soloUMIfiltering - \\
    --soloUMIdedup NoDedup \\
    --soloMultiMappers Unique \\
    --soloFeatures Gene \\
    --soloCellFilter None \\
    --soloCellReadStats Standard \\
    --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM \\
    --outSAMtype BAM Unsorted SortedByCoordinate \\
    --outSAMunmapped Within \\
    --outBAMsortingBinsN 300 \\
    """.stripIndent()
    )

    if(params.soloType == "CB_UMI_Complex"){
    scriptString.push(
    """\
    --soloType ${params.soloType} \\
    --soloCBposition ${params.soloCBposition} \\
    --soloUMIposition ${params.soloUMIposition} \\
    --soloAdapterSequence ${params.soloAdapterSequence} \\
    --soloAdapterMismatchesNmax ${params.soloAdapterMismatchesNmax}
    """.stripIndent()
    )
    }else if(params.soloType == "CB_UMI_Simple"){
    scriptString.push(
    """\
    --soloType ${params.soloType} \\
    --soloCBstart ${params.soloCBstart} \\
    --soloCBlen ${params.soloCBlen} \\
    --soloUMIstart ${params.soloUMIstart} \\
    --soloUMIlen ${params.soloUMIlen}
    """.stripIndent()
    )
    }else{
        exit 1, "soloType only support CB_UMI_Simple and CB_UMI_Complex for now."
    }

    scriptString.reverse().join()
}

process EXTRACT_BARCODE {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID), path(bam)
    tuple val(sampleID), path(merged_read1)
    tuple val(sampleID), path(merged_read2)

    output:
    tuple val(sampleID), path("*cutadapt_input.R1.fq.gz"), emit: read1
    tuple val(sampleID), path("*cutadapt_input.R2.fq.gz"), emit: read2

    script:
    def prefix   = "${sampleID}"
    def pigzThreads = Math.min(6, task.cpus)

    """
    ## Extract barcodes to a separate fastq
    reads_CB=\$(mktemp -p ./)
    temp_CB=\$(mktemp -p ./)
    tmpSAM=\$(mktemp -p ./)
    reads_CB_unique=\$(mktemp -p ./)
    
    ## covert bam to sam text with -@ parameter
    samtools view -@ ${task.cpus} $bam > \$tmpSAM
    awk '{CB=\$0; gsub(/.*CB:Z:/, "", CB); gsub(/\\t.*\$/, "", CB); if(CB!="-"){print \$1"\\t"CB}}' \$tmpSAM > \$temp_CB
    rm \$tmpSAM
    sort -u --parallel ${task.cpus} -T ./ \$temp_CB > \$reads_CB_unique
    sort -k1,1 --parallel ${task.cpus} -T ./ \$reads_CB_unique > \$reads_CB
    rm \$temp_CB \$reads_CB_unique

    tmpR1=\$(mktemp -p ./)
    tmpR1_sorted=\$(mktemp -p ./)
    tmpR1_renamed=\$(mktemp -p ./)
    zcat $merged_read1 |
        paste - - - - > \$tmpR1
    sort -t \$'\\t' -k 1,1 --parallel ${task.cpus} -T ./ \$tmpR1 > \$tmpR1_sorted
    paste \$tmpR1_sorted \$reads_CB |
        awk -F"\\t" '{print "@"\$6":"substr(\$1, 2)"\\t"substr(\$2, 47)"\\t"\$3"\\t"substr(\$4, 47)}' > \$tmpR1_renamed
    rm \$tmpR1 \$tmpR1_sorted
    
    tmpR2=\$(mktemp -p ./)
    tmpR2_sorted=\$(mktemp -p ./)
    tmpR2_renamed=\$(mktemp -p ./)
    zcat $merged_read2 |
        paste - - - - > \$tmpR2
    sort -t \$'\\t' -k 1,1 --parallel ${task.cpus} -T ./ \$tmpR1 > \$tmpR2_sorted
    paste \$tmpR2_sorted \$reads_CB |
        awk -F"\\t" '{print "@"\$6":"substr(\$1, 2)"\\t"\$2"\\t"\$3"\\t"\$4}' > \$tmpR2_renamed
    rm \$tmpR2 \$tmpR2_sorted
    rm \$reads_CB
    
    ## shuffle fastq
    tmpShuffled=\$(mktemp -p ./)
    paste \$tmpR1_renamed \$tmpR2_renamed | shuf > \$tmpShuffled
    rm \$tmpR1_renamed \$tmpR2_renamed
    
    awk -F "\\t" '{print \$1"\\n"\$2"\\n"\$3"\\n"\$4}' \$tmpShuffled |
        pigz -p $pigzThreads > ${prefix}.cutadapt_input.R1.fq.gz
    awk -F "\\t" '{print \$5"\\n"\$6"\\n"\$7"\\n"\$8}' \$tmpShuffled |
        pigz -p $pigzThreads > ${prefix}.cutadapt_input.R2.fq.gz
    rm \$tmpShuffled
    """
}

process TRIM_FASTQ {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/cutqc",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean,
        saveAs: { filename ->
        if(filename=~/readReport.json/){
            return filename
        }else{
            return null
        }
    }

    input:
    tuple val(sampleID), path(read1)
    tuple val(sampleID), path(read2)
    tuple val(sampleID), path(summary)

    output:
    tuple val(sampleID),   path("*R1.trimmed.fq.gz"), emit: read1
    tuple val(sampleID),   path("*R2.trimmed.fq.gz"), emit: read2
    tuple val(sampleID),   path("${sampleID}.readReport.json"), emit: readReport

    script:
    def prefix   = "${sampleID}"
    """
    ## cutadapt QC and trim ME adapter
    cutadapt -j $task.cpus \\
             -q 30 \\
             -m $params.trim_mLen \\
             $params.trimOpt \\
             -o ${prefix}.R1.trimmed.fq.gz -p ${prefix}.R2.trimmed.fq.gz \\
             $read1 $read2

    validBarcode_ratio=\$(awk -F"," '\$1=="Reads With Valid Barcodes"{print \$2}' $summary)
    rawReadPairs=\$(awk -F"," '\$1=="Number of Reads"{print \$2}' $summary)
    Q30_read1=\$(awk -F"," '\$1=="Q30 Bases in CB+UMI"{print \$2}' $summary)
    Q30_read2=\$(awk -F"," '\$1=="Q30 Bases in RNA read"{print \$2}' $summary)
    
    jq -n \
       --arg sampleName "${sampleID}" \\
       --arg validBarcode_ratio "\$validBarcode_ratio" \\
       --arg rawReadPairs "\$rawReadPairs" \\
       --arg Q30_read1 "\$Q30_read1" \\
       --arg Q30_read2 "\$Q30_read2" \\
       '{sampleName: \$sampleName, rawReadPairs: \$rawReadPairs, validBarcode_ratio: \$validBarcode_ratio, Q30_read1: \$Q30_read1, Q30_read2: \$Q30_read2}' > ${sampleID}.readReport.json

    """
}