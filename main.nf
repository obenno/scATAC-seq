#!/usr/bin/env nextflow

// proof of concept version
// for scATAC-seq

nextflow.enable.dsl=2

// check mandatory params
if (!params.input) { exit 1, 'Input samplesheet not specified!' }
if (!params.starIndex) { exit 1, 'STAR Index not specified!' }
//if (!params.genomeFasta) { exit 1, 'Genome Fasta not specified!' }
if (!params.genomeGTF) { exit 1, 'Genome GTF not specified!' }
if (!params.bowtie2Index) { exit 1, 'Bowtie2 Index not specified!' }
//if (!params.chromapIndex) { exit 1, 'Chromap Index not specified!' }
if (!params.whitelist) { exit 1, 'White not specified!' }

// include all processes
include { CHECK_BARCODE; EXTRACT_BARCODE; TRIM_FASTQ } from './cat_trim_fastq' addParams(params)
//include { CHROMAP } from './chromap_mapping'
include { BOWTIE2 } from './bowtie2_mapping' addParams(params)
include { LIBRARY_COMPLEXITY } from './library_complexity' addParams(params)
include { DEDUP } from './dedup' addParams(params)
include { CHECK_SATURATION } from './sequencing_saturation' addParams(params)
include { GENERATE_FRAGMENTS } from './generate_fragment' addParams(params)
include { REPORT } from './report' addParams(params)

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
    ch_starIndex = file(params.starIndex)
    //ch_chromapIndex = file(params.chromapIndex)
    ch_bowtie2Index = file(params.bowtie2Index)
    ch_genomeGTF = file(params.genomeGTF)
    ch_whitelist = Channel.fromPath(params.whitelist.split(" ").toList())

    // perform trimming and qc
    //CAT_TRIM_FASTQ(
    //    ch_fastq,
    //    ch_whitelist
    //)
    // perform bowtie2 mapping
    // bowtie genome index_base should be the same as the input fasta name
    // when performing bowtie2-build, genomeIndex needs to be set as
    // genomeIndex = "genome.fa*" in config
    //ch_blackList = file(params.blackList)
    CHECK_BARCODE(
        ch_fastq,
        ch_genomeGTF,
        ch_starIndex,
        ch_whitelist.toList()
    )

    EXTRACT_BARCODE(
        CHECK_BARCODE.out.bam,
        CHECK_BARCODE.out.merged_read1,
        CHECK_BARCODE.out.merged_read2
    )

    TRIM_FASTQ(
        EXTRACT_BARCODE.out.read1,
        EXTRACT_BARCODE.out.read2,
        CHECK_BARCODE.out.summary
    )
    
    BOWTIE2(
        TRIM_FASTQ.out.read1,
        TRIM_FASTQ.out.read2,
        ch_bowtie2Index
    )
    LIBRARY_COMPLEXITY(
        BOWTIE2.out.bam,
        BOWTIE2.out.bai
    )
    DEDUP(
        BOWTIE2.out.bam,
        BOWTIE2.out.bai
    )
    CHECK_SATURATION(
        BOWTIE2.out.bam,
        BOWTIE2.out.bai
    )    
    GENERATE_FRAGMENTS(
        BOWTIE2.out.bam,
        BOWTIE2.out.bai
    )
    REPORT(
        CAT_TRIM_FASTQ.out.readReport,
        BOWTIE2.out.mappingReport,
        DEDUP.out.dupRatio,
        LIBRARY_COMPLEXITY.out.libraryComplexity,
        CHECK_SATURATION.out.outJSON,
        GENERATE_FRAGMENTS.out.fragmentFile,
        GENERATE_FRAGMENTS.out.fragmentIndex,
        ch_genomeGTF
    )
}
