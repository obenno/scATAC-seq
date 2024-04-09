process STATS {
    tag "${sampleID}"
    label 'process_high'

    input:
    tuple val(sampleID),
          path(barcode_stats),
          path(cutadapt_json),
          path(bam_stats),
          path(dedupMetrics),
          path(inputBam),
          path(inputBamBai),
          path(fragmentFile),
          path(fragmentIndex),
          path(cells_tsv)

    output:
    tuple val(sampleID), path("${sampleID}_stats.tsv"), emit: tsv
    tuple val(sampleID), path("${sampleID}.fragments.sorted.bed.gz"), emit: fragments
    tuple val(sampleID), path("${sampleID}.fragments.sorted.bed.gz.tbi"), emit: fragmentIndex

    shell:
    '''
    total_input_readPairs=$(awk '$1=="nbr_reads:"{print $2}' !{barcode_stats})
    valid_barcode_readPairs=$(awk '$1~/nbr_reads_with_bc.*_correct_or_correctable/{print $2}' !{barcode_stats})
    with_me_q30_readPairs=$(jq -r ".read_counts.output" !{cutadapt_json})
    total_mapped_readPairs=$(grep ^SN !{bam_stats} | cut -f 2- | awk -F"\t" '$1=="reads mapped and paired:"{print $2/2}')
    duplicated_readPairs=$(awk '$0~/^DUPLICATE PAIR:/{print $NF/2}' !{dedupMetrics})
    ##uniquely_mapped_readPairs=$(samtools view -@ !{task.cpus} -F 0x4 -F 0x100 -F 0x800 -e '! [XA] && ! [SA]' !{inputBam})
    q30_readPairs=$(samtools view -c -@ !{task.cpus} -f 0x1 -f 0x2 -f 0x40 -F 0x4 -F 0x8 -F 0x100 -F 0x800 -q 30 !{inputBam})
    barcode_q30_ratio=$(awk '$1=="valid_barcode_reads_barcode_q30_ratio"{print $2}' !{barcode_stats})
    r1_q30_ratio=$(awk '$1=="valid_barcode_reads_r1_q30_ratio"{print $2}' !{barcode_stats})
    r2_q30_ratio=$(awk '$1=="valid_barcode_reads_r2_q30_ratio"{print $2}' !{barcode_stats})

    zcat !{fragmentFile} |
        awk '
        BEGIN{
            total_frag=0
            total_read=0
            mt_frag=0
            nfr_frag=0
            mono_frag=0
            m1_frag=0
            m2_frag=0
        }
        {
            total_frag+=1
            total_read+=$5
            if($1~/^chrM/){
                mt_frag+=1
            }
            if($3-$2<=147){
                nfr_frag+=1
            }else if($3-$2<=294){
                mono_frag+=1
            }

            if($5==1){
                m1_frag+=1
            }else if($5==2){
                m2_frag+=1
            }
        }
        END{
            print total_frag > "total_frag_count"
            print total_read > "total_read_count"
            print mt_frag > "mt_frag_count"
            print nfr_frag > "nfr_frag_count"
            print mono_frag > "mono_frag_count"
            print m1_frag > "read1_frag_count"
            print m2_frag > "read2_frag_count"
        }
        '

    ## pbc1 = m1_frag/total_frag
    pbc1=$(awk 'NR==FNR{m1=$1}NR>FNR{total=$1}END{print m1/total}' read1_frag_count total_frag_count)
    ## pbc1 = m1_frag/m2_frag
    pbc2=$(awk 'NR==FNR{m1=$1}NR>FNR{m2=$1}END{print m1/m2}' read1_frag_count read2_frag_count)
    ## NRF = total_frag/total_read
    NRF=$(awk 'NR==FNR{total_frag=$1}NR>FNR{total_read=$1}END{print total_frag/total_read}' total_frag_count total_read_count)

    total_fragments=$(cat total_frag_count)
    mito_fragments=$(cat mt_frag_count)
    nfr_fragments=$(cat nfr_frag_count)
    mono_fragments=$(cat mono_frag_count)

    nCells=$(wc -l !{cells_tsv} | awk '{print $1}')

    ## summarize
    echo -e "sampleID\\t!{sampleID}" > !{sampleID}_stats.tsv
    echo -e "cells\\t$nCells" >> !{sampleID}_stats.tsv
    echo -e "total_input_reads\\t$total_input_readPairs" >> !{sampleID}_stats.tsv
    echo -e "valid_barcode_reads\\t$valid_barcode_readPairs" >> !{sampleID}_stats.tsv
    echo -e "me_containing_reads\\t$with_me_q30_readPairs" >> !{sampleID}_stats.tsv
    echo -e "total_mapped_reads\\t$total_mapped_readPairs" >> !{sampleID}_stats.tsv
    ##echo -e "uniquely_mapped_reads\\t$uniquely_mapped_readPairs" >> !{sampleID}_stats.tsv
    echo -e "q30_mapped_reads\\t$q30_readPairs" >> !{sampleID}_stats.tsv
    echo -e "duplicated_reads\\t$duplicated_readPairs" >> !{sampleID}_stats.tsv
    echo -e "barcode_q30_ratio\\t$barcode_q30_ratio" >> !{sampleID}_stats.tsv
    echo -e "r1_q30_ratio\\t$r1_q30_ratio" >> !{sampleID}_stats.tsv
    echo -e "r2_q30_ratio\\t$r2_q30_ratio" >> !{sampleID}_stats.tsv
    echo -e "total_frags\\t$total_fragments" >> !{sampleID}_stats.tsv
    echo -e "mito_frags\\t$mito_fragments" >> !{sampleID}_stats.tsv
    echo -e "nfr_frags\\t$nfr_fragments" >> !{sampleID}_stats.tsv
    echo -e "mono_frags\\t$mono_fragments" >> !{sampleID}_stats.tsv
    echo -e "pbc1\\t$pbc1" >> !{sampleID}_stats.tsv
    echo -e "pbc2\\t$pbc2" >> !{sampleID}_stats.tsv
    echo -e "nrf\\t$NRF" >> !{sampleID}_stats.tsv
    '''
}