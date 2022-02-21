#!/usr/bin/env nextflow

// proof of concept version
// for scATAC-seq

nextflow.enable.dsl=2

// check mandatory params
if (!params.input) { exit 1, 'Input samplesheet not specified!' }
if (!params.genomeIndex) { exit 1, 'Genome Index not specified!' }

// include all processes
include { CAT_TRIM_FASTQ } from './cat_trim_fastq' addParams(params)
include { BOWTIE2 } from './bowtie2_mapping' addParams(params)
include { GENERATE_FRAGMENTS } from './generate_fragment' addParams(params)
include { REPORT } from './report' addParams(params)

// Reads sample list in
// code modified from nf-core/rna-seq pipeline
// https://github.com/nf-core/rnaseq
def create_fastq_channel(LinkedHashMap row) {
    def sample = [:]
    sample.name = row.sample
    sample.id = row.sample + "_" + row.rep

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }

    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }

    array = [ sample.id, sample.name, file(row.fastq_1), file(row.fastq_2) ]

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
    // perform trimming and qc
    CAT_TRIM_FASTQ( ch_fastq )
    // perform bowtie2 mapping
    // bowtie genome index_base should be the same as the input fasta name
    // when performing bowtie2-build, genomeIndex needs to be set as
    // genomeIndex = "genome.fa*" in config
    ch_genomeIndex = file(params.genomeIndex)
    ch_genomeGTF = file(params.genomeGTF)
    //ch_blackList = file(params.blackList)
    BOWTIE2(
        CAT_TRIM_FASTQ.out.read1,
        CAT_TRIM_FASTQ.out.read2,
        ch_genomeIndex
    )
    GENERATE_FRAGMENTS(
        BOWTIE2.out.bam,
        BOWTIE2.out.bai
    )
    REPORT(
        //CAT_TRIM_FASTQ.out.report,
        BOWTIE2.out.mappingReport,
        GENERATE_FRAGMENTS.out.fragmentFile,
        GENERATE_FRAGMENTS.out.fragmentIndex
    )
}
