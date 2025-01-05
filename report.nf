process REPORT {
    tag "${sample.id}"
    label 'process_high'
    publishDir "${params.outdir}/${sample.id}/final",
        mode: "${params.publish_dir_mode}",
        enabled: params.outdir as boolean
    input:
    tuple val(sample),
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
    tuple val(sample), path("${sample.id}_scATAC_report.html"), emit: report
    tuple val(sample), path("${sample.id}_combined_stats.tsv.gz"), emit: stats
    tuple val(sample), path("${sample.id}.metrics.json"), emit: metrics

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
            sampleName = "${sample.id}",
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
        output_file = "${sample.id}_scATAC_report.html"
    )
    """
}
