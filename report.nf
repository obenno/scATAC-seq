process REPORT {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/report/${sampleID}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    input:
    tuple val(sampleID), val(sampleName), path(readReport)
    tuple val(sampleID), val(sampleName), path(mappingReport)
    tuple val(sampleID), val(sampleName), path(fragmentFile)
    tuple val(sampleID), val(sampleName), path(fragmentIndex)

    output:
    tuple val(sampleID), val(sampleName), path("${sampleID}_scATAC_report.html"), emit: report

    script:
    if(params.genomeSize){
        """
        Rscript -e 'rmarkdown::render("$baseDir/bin/scATAC_report.Rmd", params = list(sampleName = "${sampleID}", fragmentFile = "${fragmentFile}", cpus = $task.cpus, refGenome = "$params.refGenome", genomeSize = "$params.genomeSize", nCount_min = "$params.nCountMin", nCount_max = "$params.nCountMax", FRiP_cutoff= "$params.FRiP_cutoff", blacklist_fraction_cutoff="$params.blacklist_fraction_cutoff", nucleosome_signal_cutoff = "$params.nucleosome_signal_cutoff", TSS_enrichment_cutoff = "$params.TSS_enrichment_cutoff", readReport = "${readReport}", bowtie2_report = "${mappingReport}"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${sampleID}_scATAC_report.html")'
        """
    }else{
        """
        Rscript -e 'rmarkdown::render("$baseDir/bin/scATAC_report.Rmd", params = list(sampleName = "${sampleID}", fragmentFile = "${fragmentFile}", cpus = $task.cpus, refGenome = "$params.refGenome", nCount_min = "$params.nCountMin", nCount_max = "$params.nCountMax", FRiP_cutoff= "$params.FRiP_cutoff", blacklist_fraction_cutoff="$params.blacklist_fraction_cutoff", nucleosome_signal_cutoff = "$params.nucleosome_signal_cutoff", TSS_enrichment_cutoff = "$params.TSS_enrichment_cutoff", readReport = "${readReport}", bowtie2_report = "${mappingReport}"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${sampleID}_scATAC_report.html")'
        """
    }
}
