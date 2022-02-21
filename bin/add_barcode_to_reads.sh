#! /usr/bin/env bash

set -euo pipefail

## This script is to simulate sinto barcode command
## the original command is too slow

bcLen=$1
## input files need to be .fq.gz
bcRead=$2
read1=$3
read2=$4

outRead1=${read1%%.fq.gz}".barcoded.fq"
outRead2=${read2%%.fq.gz}".barcoded.fq"

awk '
  FILENAME==ARGV[1]&&FNR%4==1{
    readID=substr($1,2)
    getline;
    bc[readID]=substr($1,1,'$bcLen')
  }
  FILENAME==ARGV[2]&&FNR%4==1{
    readID=substr($1,2);
    print "@"bc[readID]":"substr($1,2)" "$2 > "'$outRead1'";
    getline;
    print > "'$outRead1'";
    getline;
    print > "'$outRead1'";
    getline;
    print > "'$outRead1'";
  }
  FILENAME==ARGV[3]&&FNR%4==1{
    readID=substr($1,2);
    print "@"bc[readID]":"substr($1,2)" "$2 > "'$outRead2'";
    getline;
    print > "'$outRead2'";
    getline;
    print > "'$outRead2'";
    getline;
    print > "'$outRead2'";
  }
  ' <(zcat $bcRead) <(zcat $read1) <(zcat $read2)

pigz -p 4 $outRead1 $outRead2
