#! /usr/bin/env bash

set -euo pipefail

whitelist=$1
inputRead1=$2
bcReadLen=$3
bcLen=$4

bcOut=${inputRead1%%.fq.gz}".bc.fq"
cDNAout=${inputRead1%%.fq.gz}".cDNA_R1.fq"

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

zcat $inputRead1 |
    awk -v bcOut="$bcOut" -v cDNAout="$cDNAout" -v bcLen="$bcLen" -v bcReadLen="$bcReadLen" '
    ARGIND==1{
        wl[$1]
    }
    ARGIND==2&&FNR%4==1{
        readID=$1
        readinfo=$2
        sub("1:N", "0:N", readinfo)
        bc_line1=readID" "readinfo
        cDNA_line1=$0
        getline;
        bc=substr($1, 1, bcLen)
        bc_line2=substr($1, 1, bcReadLen)
        cDNA_line2=substr($1, bcReadLen+1)
        getline;
        line3=$0;
        getline;
        bc_line4=substr($1, 1, bcReadLen)
        cDNA_line4=substr($1, bcReadLen+1)
        if(bc in wl){
            print bc_line1 > bcOut
            print bc_line2 > bcOut
            print line3 > bcOut
            print bc_line4 > bcOut
            print cDNA_line1 > cDNAout
            print cDNA_line2 > cDNAout
            print line3 > cDNAout
            print cDNA_line4 > cDNAout
        }
        print substr($1,1, bcLen) > bcOut
        print substr($1, bcLen+1) > cDNAout
    }
    ' $whitelist_path -

pigz -p 4 $bcOut $cDNAout
