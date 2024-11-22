process REPORT {
    tag "${sampleID}"
    label 'process_high'
    publishDir "${params.outdir}/${sampleID}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    input:
    tuple val(sampleID),
          path(stats_tsv),
          path(raw_cells),
          path(raw_meta),
          path(filtered_obj),
          path(fragments),
          path(fragmentIndex),
          path(macs_peaks),
          path(saturation_json)
    path(version_json)

    output:
    tuple val(sampleID), path("${sampleID}_scATAC_report.html"), emit: report
    tuple val(sampleID), path("${sampleID}_combined_stats.tsv.gz"), emit: stats
    tuple val(sampleID), path("${sampleID}.metrics.json"), emit: metrics

    script:
    def nMem = "${task.memory.toBytes()}"
    //hgmm disabled for now
    //def rmdScript = params.refGenome == "hg38-mm10" ? "scATAC_report.hgmm.Rmd" : "scATAC_report.Rmd"
    def rmdScript = "scATAC_report.Rmd"
    """
    #! /usr/bin/env Rscript
    
    rmarkdown::render(
        "$baseDir/bin/${rmdScript}",
        params = list(
            sampleName = "${sampleID}",
            nCPUs = $task.cpus,
            nMem = "${nMem}",
            stats_tsv = "$stats_tsv",
            raw_cells = "$raw_cells",
            raw_meta = "$raw_meta",
            obj_filtered = "$filtered_obj",
            macs_peaks = "$macs_peaks",
            saturation_json = "$saturation_json",
            version_json = "${version_json}"
        ),
        intermediates_dir = getwd(),
        knit_root_dir = getwd(),
        output_dir = getwd(),
        output_file = "${sampleID}_scATAC_report.html"
    )
    """
}
