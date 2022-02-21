#! /usr/bin/env bash

set -euo pipefail

inputRead1=$1
bcLen=$2

bcOut=${inputRead1%%.fq.gz}".bc.fq"
cDNAout=${inputRead1%%.fq.gz}".cDNA_R1.fq"

zcat $inputRead1 |
    awk '
    NR%4==1{
        readID=$1
        readinfo=$2
        sub("1:N", "0:N", readinfo)
        print readID" "readinfo > "'$bcOut'";
        print $0 > "'$cDNAout'"
        getline;
        print substr($1,1,'$bcLen') > "'$bcOut'"
        print substr($1, '$bcLen'+1) > "'$cDNAout'"
        getline;
        print > "'$bcOut'"
        print > "'$cDNAout'"
        getline;
        print substr($1,1,'$bcLen') > "'$bcOut'"
        print substr($1, '$bcLen'+1) > "'$cDNAout'"
    }
    '

pigz -p 4 $bcOut $cDNAout
