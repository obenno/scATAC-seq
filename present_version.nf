process GET_VERSIONS {
    tag "Get_Version"
    label 'process_low'

    output:
    path("versions.json"), emit: json

    script:
    """
    ## fastqc version
    fastqc_version=\$(fastqc --version | awk '{print \$2}')
    ## cutadapt version
    cutadapt_version=\$(cutadapt --version | sed '/^\$/d')
    ## STAR version
    bwa_version=\$(bwa-mem2 version)
    ## samtools version
    samtools_version=\$(samtools --version | head -1 |awk '{print \$2}')
    ## bedtools version
    bedtools_version=\$(bedtools --version | awk '{print \$2}')
    cat<<-EOF > versions.json
	{
	"pipeline_version": "$workflow.manifest.version",
	"referenceGenome": "${params.refGenome}",
	"bwa_index": "${params.bwaIndex}",
	"referenceGTF": "${params.genomeGTF}",
	"bwa_version": "\$bwa_version",
	"samtools_version": "\$samtools_version",
	"bedtools_version": "\$bedtools_version"
	}
	EOF
    """
}