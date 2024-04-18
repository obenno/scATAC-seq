#! /usr/bin/env bash

set -eo pipefail

usage(){
    cat<<-EOF
Calculate saturation data from fragments file:

$(basename $0) --frags fragments.bed.gz --cells cells.tsv --peaks macs2_peaks --out out.json
EOF
}

if [[ -z $1 ]]
then
    usage
    exit 0
fi

## Read arguments
while [[ ! -z $1 ]]
do
    case $1 in
        "" | "-h" | "-help" | "--help")
            usage
            exit 0
            ;;
        "--frags")
            frags=$2
            shift 2
            ;;
        "--cells")
            cells=$2
            shift 2
            ;;
        "--peaks")
            peaks=$2
            shift 2
            ;;
        "--out")
            outputJSON=$2
            shift 2
            ;;
        "--threads")
            threads=$2
            shift 2
            ;;
    esac
done

if [[ -z $threads ]]
then
    threads=4
fi

## Expand fragment file to "reads" bed
expand_frags(){
    local frag=$1
    local readBED=$2
    zcat $frag |
        awk '
        {
            for(i=1;i<=$5;i++){
                print $1"\t"$2"\t"$3"\t"$4"-"i
            }
        }
        ' > $readBED
}

downsample_reads(){
    local readBED=$1
    awk '
    BEGIN{
        srand()
        percentage[1]=0.05
        percentage[2]=0.1
        percentage[3]=0.15
        percentage[4]=0.2
        percentage[5]=0.25
        percentage[6]=0.3
        percentage[7]=0.4
        percentage[8]=0.6
        percentage[9]=0.8
        for(i=1;i<=length(percentage);i++){
            downFile[i]="sub_"percentage[i]"_reads.bed"
        }
    }
    {
        randNumber=rand()
        for(i=1;i<=length(percentage);i++){
            if(randNumber<=percentage[i]){
                print > downFile[i]
            }
        }
    }
    ' $readBED
    ln -s $readBED sub_1_reads.bed
}

collapse_reads(){
    local readBED=$1
    local out_frag=$2
    awk '
    {
        print $1"\t"$2"\t"$3"\t"substr($4, 1, index($4, "-")-1)"\t"$4
    }
    ' $readBED |
        bedtools groupby -g 1,2,3,4 -c 5 -o count > $out_frag
}

compress_frag(){
    local frag=$1
    if [[ -z $2 ]]
    then
        local thread=1
    else
        local thread=$2
    fi

    bgzip -@ $thread $frag
    tabix -p bed $frag".gz"
}

get_metrics(){
    local cells=$1
    local peaks=$2
    local frags=$3
    get_metrics.R -p $peaks -f $frags -c $cells
}

read_bed_file=$(mktemp -p ./ XXXXXXXX.bed)
expand_frags $frags $read_bed_file
downsample_reads $read_bed_file

## collapse reads parallelly
export -f collapse_reads
parallel --link -P $threads collapse_reads ::: $(for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1}; do echo "sub_${i}_reads.bed"; done) ::: $(for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1}; do echo "sub_${i}_frags.bed"; done)

## compress frag parallelly
export -f compress_frag
parallel -P $threads compress_frag ::: $(for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1}; do echo "sub_${i}_frags.bed"; done)

for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1};
do
    fragFile="sub_${i}_frags.bed.gz"
    get_metrics $cells $peaks $fragFile |
        awk -v percentage=$i '{print percentage"\t"$0}'
done |
    jq --raw-input --slurp 'split("\n") |map(split("\t")) | .[0:-1] | map( { "percentage": .[0], "medianFragsInPeaksInCells": .[1], "medianFragsInCells": .[2], "totalFragsInCells": .[3], "total_frags": .[4], "total_reads": .[5]} )' > $outputJSON

rm $read_bed_file
rm sub_*_reads.bed sub_*_frags.bed*
