#!/bin/bash

file1=$1
file2=$2
awk -F"\t" '
    {key = $1 }
    NR == 1 {header = key}
    !(key in result) {result[key] = $0; next}
    { for (i=2; i <= NF; i++) result[key] = result[key] FS $i }
    END {
        print result[header]
        delete result[header]
        PROCINFO["sorted_in"] = "@ind_str_asc"    
        for (key in result) print result[key]
    }
' $1 $2 
