process CAT_TRIM_FASTQ {
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
    // conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    // } else {
    //     container "biocontainers/biocontainers:v1.2.0_cv1"
    // }

    input:
    tuple val(sampleID), path(read1_list), path(read2_list)
    path(gtf)
    path(starIndex)
    path(whitelist)

    output:
    //tuple val(sampleID),     path("*.chromap_aln.bam"), emit: bam
    //tuple val(sampleID),     path("*.chromap_aln.bam.bai"), emit: bai
    tuple val(sampleID),   path("*R1.trimmed.fq.gz"), emit: read1
    tuple val(sampleID),   path("*R2.trimmed.fq.gz"), emit: read2
    //tuple val(sampleID), path("*_1.merged.bc.fq.gz"), emit: read_bc
    tuple val(sampleID),   path("${sampleID}.readReport.json"), emit: readReport
    //tuple val(sampleID), path("*.read_stat.txt"), emit: report
    //tuple val(sampleID), path("*_cutqc_report.html"), emit: cutqc_report
    //path "versions.yml"                       , emit: versions

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

    scriptString.push(
    """
    ## Extract barcodes to a separate fastq
    reads_CB=\$(mktemp -p ./)
    temp_CB=\$(mktemp -p ./)
    samtools view -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam |
        awk '{CB=\$0; gsub(/.*CB:Z:/, "", CB); gsub(/\\t.*\$/, "", CB); if(CB!="-"){print \$1"\\t"CB}}' > \$temp_CB
    sort -u --parallel ${task.cpus} -T ./ \$temp_CB > \$reads_CB
    rm \$temp_CB
    ##readsID=\$(mktemp -p ./)
    ##cut -f 1 \$reads_CB > \$readsID
    ##tmp_whitelist=\$(mktemp -p ./)
    ##cut -f 2 \$reads_CB | sort -u --parallel ${task.cpus} > \$tmp_whitelist
    ## prepare input reads for chromap
    ## use the last 90 nt of R1
    zcat ${prefix}_1.merged.fq.gz |
        awk '
        ARGIND==1{
            a[\$1]=\$2
        }
        ARGIND==2 && FNR%4==1{
            if(substr(\$1, 2) in a){
                print "@"a[substr(\$1,2)]":"substr(\$0,2);
                getline;
                print substr(\$1, 90);
                getline;
                print;
                getline;
                print substr(\$1, 90)
            }
        }
        ' \$reads_CB - |
        pigz -p $pigzThreads > ${prefix}.cutadapt_input.R1.fq.gz

    zcat ${prefix}_2.merged.fq.gz |
        awk '
        ARGIND==1{
            a[\$1]=\$2
        }
        ARGIND==2 && FNR%4==1{
            if(substr(\$1, 2) in a){
                print "@"a[substr(\$1,2)]":"substr(\$0,2);
                getline;
                print;
                getline;
                print;
                getline;
                print;
            }
        }
        ' \$reads_CB - |
        pigz -p $pigzThreads > ${prefix}.cutadapt_input.R2.fq.gz

    ##seqtk subseq ${prefix}_2.merged.fq.gz \$readsID |
    ##    pigz -p $pigzThreads > ${prefix}.cutadapt_input.R2.fq.gz

    ## cutadapt QC and trim ME adapter
    cutadapt -j $task.cpus \\
             -q 30 \\
             -m 50 \\
             $params.trimOpt \\
             -o ${prefix}.R1.trimmed.fq.gz -p ${prefix}.R2.trimmed.fq.gz \\
             ${prefix}.cutadapt_input.R1.fq.gz ${prefix}.cutadapt_input.R2.fq.gz

    validBarcode_ratio=\$(awk -F"," '\$1=="Reads With Valid Barcodes"{print \$2}' ${prefix}.Solo.out/Gene/Summary.csv)
    rawReadPairs=\$(awk -F"," '\$1=="Number of Reads"{print \$2}' ${prefix}.Solo.out/Gene/Summary.csv)
    Q30_read1=\$(awk -F"," '\$1=="Q30 Bases in CB+UMI"{print \$2}' ${prefix}.Solo.out/Gene/Summary.csv)
    Q30_read2=\$(awk -F"," '\$1=="Q30 Bases in RNA read"{print \$2}' ${prefix}.Solo.out/Gene/Summary.csv)
    
    jq -n \
       --arg sampleName "${sampleID}" \\
       --arg validBarcode_ratio "\$validBarcode_ratio" \\
       --arg rawReadPairs "\$rawReadPairs" \\
       --arg Q30_read1 "\$Q30_read1" \\
       --arg Q30_read2 "\$Q30_read2" \\
       '{sampleName: \$sampleName, rawReadPairs: \$rawReadPairs, validBarcode_ratio: \$validBarcode_ratio, Q30_read1: \$Q30_read1, Q30_read2: \$Q30_read2}' > ${sampleID}.readReport.json

    """.stripIndent()
    )
    scriptString.reverse().join()
}
