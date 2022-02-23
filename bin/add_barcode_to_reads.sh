#! /usr/bin/env bash

set -euo pipefail

## This script is to simulate sinto barcode command
## the original command is too slow

## Additionally, whitelist is used
## here to filter barcodes

## metrics will be written to stdout as json
sampleName=$1
whitelist=$2
bcLen=$3
## input files need to be .fq.gz
bcRead=$4
read1=$5
read2=$6

outbc=${bcRead%%.bc.fq.gz}".whitelist_bc.fq"
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

awk -v bcLen="$bcLen" -v outbc="$outbc" -v outRead1="$outRead1" -v outRead2="$outRead2" '
  FILENAME==ARGV[1]{
    whitelist[$1]
  }
  FILENAME==ARGV[2]&&FNR%4==1{
    readID=substr($1,2)
    firstLine=$0
    getline;
    bc_seq=substr($1,1, bcLen)
    if(bc_seq in whitelist){
      bc[readID]=bc_seq
      print firstLine > outbc
      print bc_seq > outbc
      getline;
      print > outbc
      getline
      pirnt substr($1, 1, bcLen) > outbc
    }
  }
  FILENAME==ARGV[3]&&FNR%4==1{
    readID=substr($1,2);
    if(readID in bc){
      print "@"bc[readID]":"substr($1,2)" "$2 > outRead1;
      getline;
      print > outRead1;
      getline;
      print > outRead1;
      getline;
      print > outRead1;
    }
  }
  FILENAME==ARGV[4]&&FNR%4==1{
    readID=substr($1,2);
    if(readID in bc){
      print "@"bc[readID]":"substr($1,2)" "$2 > outRead2;
      getline;
      print > outRead2;
      getline;
      print > outRead2;
      getline;
      print > outRead2;
    }
  }
  ' $whitelist_path <(zcat $bcRead) <(zcat $read1) <(zcat $read2)

validBarcode_count=$(awk 'NR%4==1{split($1, tmp, ":"); barcode=substr(tmp[1],2); print barcode}' $outRead1 | sort -u | wc -l)
validBarcode_readCount=$(awk 'END{print NR/4}' $outRead1)
rawReadPairs=$(zcat $read1 | awk 'END{print NR/4*2}')
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
