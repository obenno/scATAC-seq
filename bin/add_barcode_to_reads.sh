#! /usr/bin/env bash

set -euo pipefail

## This script is to simulate sinto barcode command
## the original command is too slow

## Additionally, whitelist is used
## here to filter barcodes

## This script assumes bc sequence is
## from R1, from position 1-bcReadLen,
## and from 1-bcLen is barcode sequence
## (barcode reads contains sequences other than barcode)

## metrics will be written to stdout as json
sampleName=$1
whitelist=$2
bcReadLen=$3
bcLen=$4
## input files need to be .fq.gz
read1=$5
read2=$6
threads=$7

outbc=${read1%%.fq.gz}".whitelist_bc.fq"
outRead1=${read1%%.fq.gz}".barcoded.fq"
outRead2=${read2%%.fq.gz}".barcoded.fq"

if [[ ! -f $whitelist ]]
then
    echo "Please provide whitelist file"
    exit 1
elif [[ $whitelist =~ \.gz$ ]]
then
    whitelist_path=$(mktemp -p .)
    gunzip -c $whitelist > $whitelist_path
else
    whitelist_path=$whitelist
fi

## sort R1 and R2
read1_temp=$(mktemp -p ./)
zcat $read1 | paste - - - - |
    sort -k 1,1 --parallel=$threads -S 10% -T ./ > $read1_temp

read2_temp=$(mktemp -p ./)
zcat $read2 | paste - - - - |
    sort -k 1,1 --parallel=$threads -S 10% -T ./ > $read2_temp

paste -d "\t" $read1_temp $read2_temp |
    awk -F "\t" -v bcLen="$bcLen" -v bcReadLen="$bcReadLen" -v outbc="$outbc" -v outRead1="$outRead1" -v outRead2="$outRead2" '
    ARGIND==1{
        wl[$1]
    }
    ARGIND==2{
        split($1, tmp1, " ")
        read1_id=tmp1[1]
        split($5, tmp2, " ")
        read2_id=tmp2[1]
        if(read1_id!=read2_id){
            print "R1 and R2 are not concordant, please check fastq input files." > "/dev/stderr"
            exit 1
        }
        bc_seq=substr($2, 1, bcLen)
        bc_qual=substr($4, 1, bcLen)
        cDNA_seq=substr($2, bcReadLen+1)
        cDNA_qual=substr($4, bcReadLen+1)
        if(bc_seq in wl){
            print "@"bc_seq":"substr($1,2)"\n"bc_seq"\n"$3"\n"bc_qual > outbc
            print "@"bc_seq":"substr($1,2)"\n"cDNA_seq"\n"$3"\n"cDNA_qual > outRead1
            print "@"bc_seq":"substr($5,2)"\n"$6"\n"$7"\n"$8 > outRead2
        }
    }
    ' $whitelist_path -
rm $read1_temp $read2_temp

##awk -v bcLen="$bcLen" -v outbc="$outbc" -v outRead1="$outRead1" -v outRead2="$outRead2" '
##  FILENAME==ARGV[1]{
##    whitelist[$1]
##  }
##  FILENAME==ARGV[2]&&FNR%4==1{
##    readID=substr($1,2)
##    firstLine=$0
##    getline;
##    bc_seq=substr($1,1, bcLen)
##    if(bc_seq in whitelist){
##      bc[readID]=bc_seq
##      print firstLine > outbc
##      print bc_seq > outbc
##      getline;
##      print > outbc
##      getline
##      pirnt substr($1, 1, bcLen) > outbc
##    }
##  }
##  FILENAME==ARGV[3]&&FNR%4==1{
##    readID=substr($1,2);
##    if(readID in bc){
##      print "@"bc[readID]":"substr($1,2)" "$2 > outRead1;
##      getline;
##      print > outRead1;
##      getline;
##      print > outRead1;
##      getline;
##      print > outRead1;
##    }
##  }
##  FILENAME==ARGV[4]&&FNR%4==1{
##    readID=substr($1,2);
##    if(readID in bc){
##      print "@"bc[readID]":"substr($1,2)" "$2 > outRead2;
##      getline;
##      print > outRead2;
##      getline;
##      print > outRead2;
##      getline;
##      print > outRead2;
##    }
##  }
##  ' $whitelist_path <(zcat $bcRead) <(zcat $read1) <(zcat $read2)

validBarcode_count=$(awk 'NR%4==1{split($1, tmp, ":"); barcode=substr(tmp[1],2); print barcode}' $outRead1 | sort -u --parallel=$threads -T ./| wc -l)
validBarcode_readCount=$(awk 'END{print NR/4}' $outRead1)
rawReadPairs=$(zcat $read1 | awk 'END{print NR/4}')
Q30_bcRead=$(check_q30.awk $outbc)
Q30_read1=$(check_q30.awk $outRead1)
Q30_read2=$(check_q30.awk $outRead2)

pigz -p 4 $outbc $outRead1 $outRead2

## write json metrics
jq -n \
   --arg sampleName "$sampleName" \
   --arg validBarcode_count "$validBarcode_count" \
   --arg validBarcode_readCount "$validBarcode_readCount" \
   --arg rawReadPairs "$rawReadPairs" \
   --arg Q30_bcRead "$Q30_bcRead" \
   --arg Q30_read1 "$Q30_read1" \
   --arg Q30_read2 "$Q30_read2" \
   '{sampleName: $sampleName, rawReadPairs: $rawReadPairs, validBarcode_count: $validBarcode_count, validBarcode_readCount: $validBarcode_readCount, Q30_bcRead: $Q30_bcRead, Q30_read1: $Q30_read1, Q30_read2: $Q30_read2}'
