#!/usr/bin/env bash

#Author: Ian Rambo
#Purpose: Create a Gold OID to Sample tab-separated mapping file from IMG

img_directory=$1
gold_map=$2


printf "GOLD_OID\tSAMPLE_ID\n" > $gold_map

for config in $(find $img_directory -type f -name "*.config"); do
    paste -d'\t' <(grep 'analysis_project_id'  $config | cut -f2 -d' ') <(grep taxon_name $config | rev | cut -f1 -d' ' | rev | awk '{gsub(/^ +| +$/,"")} {print $0}') >> $gold_map;
done
