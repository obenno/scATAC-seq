process LIBRARY_COMPLEXITY {
    // count mitocondrial reads propertion in the raw alignment
    tag "${sampleID}"
    label 'process_low'

    input:
    tuple val(sampleID), path(inputBam)
    tuple val(sampleID), path(inputBamBai)

    output:
    tuple val(sampleID), path("${sampleID}.libComplexity"), emit: libraryComplexity

    shell:
    '''
    samtools view -@ !{task.cpus} -u -f 0x2 -F 0x4 -F 0x8 -F 0x100 -F 0x800 -F 0x400 -q 30 !{inputBam} |
        samtools sort -@ !{task.cpus} -n -u |
        bedtools bamtobed -bedpe -i - |
        awk 'BEGIN{OFS="\\t"}{print $1, $2, $4, $6, $9, $10}' |
        grep -v !{params.mito_chr_label} | sort --parallel=!{task.cpus} -S 10% -T ./ | uniq -c |
        awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{print "mt\\tm0\\tm1\\tm2\\tm0/mt\\tm1/m0\\tm1/m2"; if(mt==0){nrf="NA"}else{nrf=m0/mt};if(m0==0){pbc1="NA"}else{pbc1=m1/m0}; if(m2==0){pbc2="NA"}else{pbc2=m1/m2};printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,nrf,pbc1,pbc2}' > !{sampleID}.libComplexity
    '''
}
