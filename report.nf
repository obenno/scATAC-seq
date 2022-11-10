process REPORT {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/report/${sampleID}",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    input:
    tuple val(sampleID), val(sampleName), path(readReport)
    tuple val(sampleID), val(sampleName), path(mappingReport)
    tuple val(sampleID), val(sampleName), val(dupRatio)
    tuple val(sampleID), val(sampleName), path(libraryComplexity)
    tuple val(sampleID), val(sampleName), path(saturation_json)
    tuple val(sampleID), val(sampleName), path(fragmentFile)
    tuple val(sampleID), val(sampleName), path(fragmentIndex)
    path refGTF

    output:
    tuple val(sampleID), val(sampleName), path("${sampleID}_scATAC_report.html"), emit: report

    script:
    def nMem = "${task.memory.toBytes()}"
    if(params.genomeSize){
        """
        Rscript -e 'rmarkdown::render("$baseDir/bin/scATAC_report.Rmd", params = list(sampleName = "${sampleID}", fragmentFile = "${fragmentFile}", nCPUs = $task.cpus, nMem = "${nMem}", refGenome = "$params.refGenome", refGTF = "${refGTF}", genomeSize = "$params.genomeSize", nCount_min = "$params.nCountMin", nCount_max = "$params.nCountMax", FRiP_cutoff= "$params.FRiP_cutoff", blacklist_fraction_cutoff="$params.blacklist_fraction_cutoff", nucleosome_signal_cutoff = "$params.nucleosome_signal_cutoff", TSS_enrichment_cutoff = "$params.TSS_enrichment_cutoff", readReport = "${readReport}", bowtie2_report = "${mappingReport}", saturation_json = "$saturation_json", libraryComplexity = "$libraryComplexity", dupRatio = "$dupRatio"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${sampleID}_scATAC_report.html")'
        """
    }else{
        """
        Rscript -e 'rmarkdown::render("$baseDir/bin/scATAC_report.Rmd", params = list(sampleName = "${sampleID}", fragmentFile = "${fragmentFile}", nCPUs = $task.cpus, nMem = "${nMem}", refGenome = "$params.refGenome", refGTF = "${refGTF}", nCount_min = "$params.nCountMin", minCell = "$params.minCell", minFeature = "$params.minFeature", nCount_max = "$params.nCountMax", FRiP_cutoff= "$params.FRiP_cutoff", blacklist_fraction_cutoff="$params.blacklist_fraction_cutoff", nucleosome_signal_cutoff = "$params.nucleosome_signal_cutoff", TSS_enrichment_cutoff = "$params.TSS_enrichment_cutoff", readReport = "${readReport}", bowtie2_report = "${mappingReport}", saturation_json = "$saturation_json", libraryComplexity = "$libraryComplexity", dupRatio = "$dupRatio"), intermediates_dir = getwd(), knit_root_dir = getwd(), output_dir = getwd(), output_file = "${sampleID}_scATAC_report.html")'
        """
    }
}
