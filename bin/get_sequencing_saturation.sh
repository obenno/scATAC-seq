#! /usr/bin/env bash

## This script is to generate data for
## sequencing saturation curve

## The sequencing saturation inrespects cells

## Input is BOWTIE2 mapping BAM file
## dependency: samtools, gawk, sort, sinto

## No whitelist is required in this step, since the input reads
## of BOWTIE2 were already filtered according to whitelist

set -euo pipefail

##usage="Usage: $(basename $0) inputBAM thread refGenome outputJSON
##refGenome: hg38, hg19, mm10 or mm9
##"

usage="Usage: $(basename $0) inputBAM thread outputJSON"

if [[ $# -eq 0 ]] || [[ $1 == "-h" ]] || [[ $1 == "-help" ]] || [[ $1 == "--help" ]]
then
    echo "$usage"
    exit 0
fi

inputBAM=$1
thread=$2
##refGenome=$3
outputJSON=$3

##if [[ -z $refGenome ]]
##then
##    echo "refGenome not provided."
##    exit 1
##elif [[ $refGenome != "hg38" ]] && [[ $refGenome != "hg19" ]] && [[ $refGenome != "mm9" ]] && [[ $refGenome != "mm10" ]]
##then
##    echo "supported genomes: hg38, hg19, mm10, mm9"
##    exit 1
##fi

if [[ ! -f $inputBAM ]]
then
    echo "$inputBAM not found."
    exit 1
fi

function calc_saturation {
    local subsampledBAM=$(mktemp -p ./ subsampled.XXXXXXXX.bam)
    local subsampledFragment=$(mktemp -p ./ frag.XXXXXXXX.bed)
    ##local readInfo=$(mktemp -p ./)
    ## percentage is  0.0 ≤ FLOAT ≤ 1.0
    local percentage=$1

    samtools view -@ $thread -u --subsample $percentage --subsample-seed 1324 $inputBAM |
        samtools sort -@ $thread > $subsampledBAM
    samtools index $subsampledBAM
    sinto fragments -b $subsampledBAM -p $thread -f $subsampledFragment --barcode_regex "[^:]*" ##--collapse_within
    seqDep=$(awk '{n+=$5}END{print n}' $subsampledFragment)
    uniqueFrag=$(wc -l $subsampledFragment | awk '{print $1}')
    ##sort -k1,1 -k2,2n $subsampledFragment | bgzip -@ $thread -c > $subsampledFragment".gz"
    ##tabix -p bed $subsampledFragment".gz"
    ##medianFrag=$(report_cells.R -i $subsampledFragment".gz" -t $thread -g $refGenome)
    printf "%s\t%s\t%s\n" $percentage $seqDep $uniqueFrag
    rm $subsampledBAM $subsampledBAM".bai"
    rm $subsampledFragment ##$subsampledFragment".gz" $subsampledFragment".gz.tbi"
}

## Start time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation started at %s %s\n" $currentTime

for i in {0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.6,0.8,1};
do
    calc_saturation $i
done |
    jq --raw-input --slurp 'split("\n") |map(split("\t")) | .[0:-1] | map( { "percentage": .[0], "seqDep": .[1], "uniqueFrags": .[2] } )' > $outputJSON

## End time
currentTime=$(date +"%F %T")
>&2 printf "Calculating sequencing saturation ended at %s %s\n" $currentTime
