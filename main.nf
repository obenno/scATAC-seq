#!/usr/bin/env nextflow

// proof of concept version
// for scATAC-seq

nextflow.enable.dsl=2

// check mandatory params
if (!params.input) { exit 1, 'Input samplesheet not specified!' }
if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }
if (!params.bwaIndex) { exit 1, 'BWA Index not specified!' }
if (!params.whitelist) { exit 1, 'Whitelist not specified!' }

// include all processes
include { CAT_FASTQ; CHECK_BARCODE; TRIM_FASTQ } from './cat_trim_fastq'
include { CAT_FASTQ_10X; CHECK_BARCODE_10X; } from './cat_trim_fastq'
include { BWA_MAPPING } from './bwa_mapping'
include { MULTIQC } from './multiqc'
include { STATS } from './stats'
include { SIGNAC } from './signac_process'
include { DEDUP } from './dedup'
include { CHECK_SATURATION } from './sequencing_saturation'
include { GENERATE_FRAGMENTS } from './generate_fragment'
include { GET_VERSIONS } from './present_version'
include { REPORT } from './report'

// Reads sample list in
// code modified from nf-core/rna-seq pipeline
// https://github.com/nf-core/rnaseq
def create_fastq_channel_TB(LinkedHashMap row) {
    def sample = [:]

    if(!row.sample){
        exit 1, "ERROR: Please check sample list header and ensure that \"sample\" exists!\n"
    }else{
        sample.id = row.sample
    }

    if(!row.expected_cells && params.topCells){
        exit 1, "ERROR: Please check input samplesheet -> Please specify expected_cells column"
    }else if(!params.topCells){
        sample.expected_cells = 0
    }else{
        def expected_cells = row.expected_cells as Integer
        assert expected_cells.getClass() == Integer
        sample.expected_cells = expected_cells
    }

    if(!row.fastq_1){
        exit 1, "ERROR: Please check sample list header and ensure that \"fastq_1\" exists!\n"
    }
    if(!row.fastq_2){
        exit 1, "ERROR: Please check sample list header and ensure that \"fastq_2\" exists!\n"
    }
    
    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }

    //array = [ sample.id, sample.name, file(row.fastq_1), file(row.fastq_2) ]
    array = [ sample, file(row.fastq_1), file(row.fastq_2) ]

    return array
}

def create_fastq_channel_10X(LinkedHashMap row) {
    def sample = [:]

    if(!row.sample){
        exit 1, "ERROR: Please check sample list header and ensure that \"sample\" exists!\n"
    }else{
        sample.id = row.sample
    }

    if(!row.expected_cells && params.topCells){
        exit 1, "ERROR: Please check input samplesheet -> Please specify expected_cells column"
    }else if(!params.topCells){
        sample.expected_cells = 0
    }else{
        def expected_cells = row.expected_cells as Integer
        assert expected_cells.getClass() == Integer
        sample.expected_cells = expected_cells
    }

    if(!row.fastq_1){
        exit 1, "ERROR: Please check sample list header and ensure that \"fastq_1\" exists!\n"
    }
    if(!row.fastq_2){
        exit 1, "ERROR: Please check sample list header and ensure that \"fastq_2\" exists!\n"
    }
    if(!row.fastq_3){
        exit 1, "ERROR: Please check sample list header and ensure that \"fastq_3\" exists!\n"
    }

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }

    if (!file(row.fastq_3).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 3 FastQ file does not exist!\n${row.fastq_3}"
    }

    array = [ sample, file(row.fastq_1), file(row.fastq_2), file(row.fastq_3) ]

    return array
}

workflow {
    main:
    ch_bwaIndex = file(params.bwaIndex + "*")
    ch_genomeGTF = file(params.genomeGTF)

    if(params.platform == "TB"){

        // TB used three barcode files
        ch_whitelist = Channel.fromPath(params.whitelist.split(" ").toList())

        Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map{ create_fastq_channel_TB(it) }
        .groupTuple(by: 0)
        .set{ ch_fastq }

        CAT_FASTQ( ch_fastq )

        CAT_FASTQ.out.read1
        .join(CAT_FASTQ.out.read2, by:[0])
        .set{ ch_merged_fastq }


        CHECK_BARCODE(
            ch_merged_fastq,
            ch_whitelist.toList()
        )

        ch_barcode_trimmed_fastq = CHECK_BARCODE.out.read1.join(CHECK_BARCODE.out.read2, by:[0])

        TRIM_FASTQ(
            ch_barcode_trimmed_fastq
        )

        ch_me_trimmed_fastq = TRIM_FASTQ.out.read1.join(TRIM_FASTQ.out.read2, by:[0])

        BWA_MAPPING(
            ch_me_trimmed_fastq,
            ch_bwaIndex
        )
    
        ch_mapping_bam = BWA_MAPPING.out.bam.join(BWA_MAPPING.out.bai, by:[0])

        DEDUP(
            ch_mapping_bam
        )

        CAT_FASTQ.out.read1
        .join(CAT_FASTQ.out.read2, by:[0])
        .join(TRIM_FASTQ.out.read1, by:[0])
        .join(TRIM_FASTQ.out.read2, by:[0])
        .join(TRIM_FASTQ.out.report_JSON, by:[0])
        .join(DEDUP.out.bam_stats, by:[0])
        .set{ ch_read_stats }

        MULTIQC(
            ch_read_stats
        )

        GENERATE_FRAGMENTS(
            ch_mapping_bam
        )

        GENERATE_FRAGMENTS.out.fragmentFile
        .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
        .set{ ch_fragment }

        SIGNAC(
            ch_fragment,
            ch_genomeGTF
        )

        SIGNAC.out.raw_cells
        .join(SIGNAC.out.macs_peaks, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentFile, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
        .set{ ch_saturation_input }

        CHECK_SATURATION(
            ch_saturation_input
        )

        // STATS input channel
        CHECK_BARCODE.out.stats
        .join(TRIM_FASTQ.out.report_JSON, by:[0])
        .join(DEDUP.out.bam_stats, by:[0])
        .join(DEDUP.out.dedup_metrics, by:[0])
        .join(BWA_MAPPING.out.bam, by:[0])
        .join(BWA_MAPPING.out.bai, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentFile, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
        .join(SIGNAC.out.raw_cells, by:[0])
        .set{ ch_stats_input }

    }else if(params.platform == "10X"){

        ch_whitelist = file(params.whitelist)

        Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map{ create_fastq_channel_10X(it) }
        .groupTuple(by: 0)
        .set{ ch_fastq }

        CAT_FASTQ_10X( ch_fastq )

        CAT_FASTQ_10X.out.read2
        .join(CAT_FASTQ_10X.out.read1, by:[0])
        .join(CAT_FASTQ_10X.out.read3, by:[0])
        .set{ ch_merged_fastq }

        CHECK_BARCODE_10X(
            ch_merged_fastq,
            ch_whitelist
        )

        ch_barcode_trimmed_fastq = CHECK_BARCODE_10X.out.read1.join(CHECK_BARCODE_10X.out.read2, by:[0])

        TRIM_FASTQ(
            ch_barcode_trimmed_fastq
        )

        ch_me_trimmed_fastq = TRIM_FASTQ.out.read1.join(TRIM_FASTQ.out.read2, by:[0])

        BWA_MAPPING(
            ch_me_trimmed_fastq,
            ch_bwaIndex
        )
    
        ch_mapping_bam = BWA_MAPPING.out.bam.join(BWA_MAPPING.out.bai, by:[0])

        DEDUP(
            ch_mapping_bam
        )

        CAT_FASTQ_10X.out.read1
        .join(CAT_FASTQ_10X.out.read3, by:[0])
        .join(TRIM_FASTQ.out.read1, by:[0])
        .join(TRIM_FASTQ.out.read2, by:[0])
        .join(TRIM_FASTQ.out.report_JSON, by:[0])
        .join(DEDUP.out.bam_stats, by:[0])
        .set{ ch_read_stats }

        MULTIQC(
            ch_read_stats
        )

        GENERATE_FRAGMENTS(
            ch_mapping_bam
        )

        GENERATE_FRAGMENTS.out.fragmentFile
        .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
        .set{ ch_fragment }

        SIGNAC(
            ch_fragment,
            ch_genomeGTF
        )

        SIGNAC.out.raw_cells
        .join(SIGNAC.out.macs_peaks, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentFile, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
        .set{ ch_saturation_input }

        CHECK_SATURATION(
            ch_saturation_input
        )
        
        // STATS input channel
        CHECK_BARCODE_10X.out.stats
        .join(TRIM_FASTQ.out.report_JSON, by:[0])
        .join(DEDUP.out.bam_stats, by:[0])
        .join(DEDUP.out.dedup_metrics, by:[0])
        .join(BWA_MAPPING.out.bam, by:[0])
        .join(BWA_MAPPING.out.bai, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentFile, by:[0])
        .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
        .join(SIGNAC.out.raw_cells, by:[0])
        .set{ ch_stats_input }
        
    }

    STATS(
        ch_stats_input
    )

    GET_VERSIONS()
    
    // join report input
    STATS.out.tsv
    .join(SIGNAC.out.raw_cells, by:[0])
    .join(SIGNAC.out.raw_meta, by:[0])
    .join(SIGNAC.out.obj, by:[0])
    .join(GENERATE_FRAGMENTS.out.fragmentFile, by:[0])
    .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
    .join(SIGNAC.out.macs_peaks, by:[0])
    .join(CHECK_SATURATION.out.outJSON, by:[0])
    .set{ ch_report_input }

    REPORT(
        ch_report_input,
        GET_VERSIONS.out.json
    )
}
