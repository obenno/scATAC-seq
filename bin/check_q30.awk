#! /usr/bin/env -S awk -f

## The code is from gnu awk's online manual
## https://www.gnu.org/software/gawk/manual/html_node/Ordinal-Functions.html


BEGIN {
    _ord_init()
    q30=0
    n=0
}

function _ord_init(    low, high, i, t){
    low = sprintf("%c", 7) # BEL is ascii 7
    if (low == "\a") {    # regular ascii
        low = 0
        high = 127
    } else if (sprintf("%c", 128 + 7) == "\a") {
        # ascii, mark parity
        low = 128
            high = 255
    } else {        # ebcdic(!)
        low = 0
        high = 255
    }

    for (i = low; i <= high; i++) {
        t = sprintf("%c", i)
        _ord_[t] = i
    }
}

function ord(str,    c){
    # only first character is of interest
    c = substr(str, 1, 1)
    return _ord_[c]
}

function chr(c){
    # force c to be numeric by adding 0
    return sprintf("%c", c + 0)
}

NR%4==1{
    getline;
    getline;
    getline;
    for(i=1;i<=length($1);i++){
        n++
        if(ord($1)-33>=30){
            q30++
        }
    }
}

END{
    printf "%.2f%\n", q30/n*100
}
