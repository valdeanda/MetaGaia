#!/bin/bash

cd ../../data/IMG_data/

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
' $1 $2  > $1.mapping.tsv
sed  -e "s/\r//g"  $1.mapping.tsv |  sed 's/ /\t/g' >  $1.mapping_sorted.tsv

cut -f 1,2 $1.mapping_sorted.tsv  | awk '{if (!$2) {print $1, "NoBin"} else {print $1, $2}}' > mappingFile1.tab 
#again, remove matching characters and replace by tabs

sed -i -e "s/\r//g" mappingFile1.tab 
sed -i 's/ /\t/g' mappingFile1.tab