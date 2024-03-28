#!/usr/bin/env nextflow

// proof of concept version
// for scATAC-seq

nextflow.enable.dsl=2

// check mandatory params
if (!params.input) { exit 1, 'Input samplesheet not specified!' }
if (!params.bwaIndex) { exit 1, 'STAR Index not specified!' }
//if (!params.genomeFasta) { exit 1, 'Genome Fasta not specified!' }
if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }
if (!params.bwaIndex) { exit 1, 'BWA Index not specified!' }
if (!params.whitelist) { exit 1, 'White not specified!' }

// include all processes
include { CHECK_BARCODE; TRIM_FASTQ } from './cat_trim_fastq'
include { BWA_MAPPING } from './bwa_mapping'
include { MULTIQC } from './multiqc'
include { STATS } from './stats'
include { SIGNAC } from './signac_process'
include { DEDUP } from './dedup'
include { CHECK_SATURATION } from './sequencing_saturation'
include { GENERATE_FRAGMENTS } from './generate_fragment'
include { REPORT } from './report'

// Reads sample list in
// code modified from nf-core/rna-seq pipeline
// https://github.com/nf-core/rnaseq
def create_fastq_channel(LinkedHashMap row) {
    def sample = [:]
    //sample.name = row.sample
    //sample.id = row.sample + "_" + row.rep
    sample.id = row.sample
    
    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }

    //array = [ sample.id, sample.name, file(row.fastq_1), file(row.fastq_2) ]
    array = [ sample.id, file(row.fastq_1), file(row.fastq_2) ]

    return array
}

workflow {
    main:
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map{ create_fastq_channel(it) }
        .groupTuple(by: 0)
        .set{ ch_fastq }
    ch_bwaIndex = file(params.bwaIndex + "*")
    ch_genomeGTF = file(params.genomeGTF)
    ch_whitelist = Channel.fromPath(params.whitelist.split(" ").toList())

    CHECK_BARCODE(
        ch_fastq,
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
    
    ch_mapping_bam = BWA_MAPPING.out.bam.join(BWA_MAPPING.out.bai)
    DEDUP(
        ch_mapping_bam
    )
    
    CHECK_BARCODE.out.merged_r1
    .join(CHECK_BARCODE.out.merged_r2, by:[0])
    .join(TRIM_FASTQ.out.read1, by:[0])
    .join(TRIM_FASTQ.out.read2, by:[0])
    .join(TRIM_FASTQ.out.report_JSON, by:[0])
    .join(DEDUP.out.bam_stats, by:[0])
    .set{ ch_read_stats }

    MULTIQC(
        ch_read_stats
    )

    CHECK_SATURATION(
        ch_mapping_bam
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

    // STATS input channel
    CHECK_BARCODE.out.stats
    .join(CHECK_BARCODE.out.read1, by:[0])
    .join(CHECK_BARCODE.out.read2, by:[0])
    .join(TRIM_FASTQ.out.report_JSON, by:[0])
    .join(DEDUP.out.bam_stats, by:[0])
    .join(DEDUP.out.dedup_metrics, by:[0])
    .join(BWA_MAPPING.out.bam, by:[0])
    .join(BWA_MAPPING.out.bai, by:[0])
    .join(GENERATE_FRAGMENTS.out.fragmentFile, by:[0])
    .join(GENERATE_FRAGMENTS.out.fragmentIndex, by:[0])
    .join(SIGNAC.out.raw_cells, by:[0])
    .set{ ch_stats_input }

    STATS(
        ch_stats_input
    )
    
    // join report input
    STATS.out.tsv
    .join(SIGNAC.out.raw_cells, by:[0])
    .join(SIGNAC.out.raw_meta, by:[0])
    .join(SIGNAC.out.obj, by:[0])
    .join(SIGNAC.out.macs_peaks, by:[0])
    .join(CHECK_SATURATION.out.outJSON, by:[0])
    .set{ ch_report_input }

    REPORT(
        ch_report_input
    )
}
