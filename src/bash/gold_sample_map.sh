#!/usr/bin/env bash

#Author: Ian Rambo
#Created: September 12, 2021
#Purpose: Create a Gold OID to Sample tab-separated mapping file from IMG

usage="
$(basename "$0"): Create a Gold OID to Sampling_Site tab-separated mapping file from IMG

where:

   -h --- show this help message
   -i --- directory containing IMG analysis files
   -o --- output gold_oid --> sampling_site map file
   -e --- additional string to append to the Sampling_Site

The script will search for .config files using find within the directory specified with -i

If you want to search multiple directories, pass the base directory containing
multiple downloaded IMG analysis directories.
    "

while getopts ':hi:o:e:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    i) img_directory=${OPTARG};;
    o) gold_map=${OPTARG};;
    e) extension=${OPTARG};;

    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))
#=============================================================================
test -d $img_directory || { echo "ERROR: IMG directory $img_directory does not exist, exiting..."; exit 1; }

printf "GOLD_OID\tSampling_Site\n" > $gold_map

for config in $(find $img_directory -type f -name "*.config"); do

    if [[ -s $config ]]; then

        echo "parsing: $config"

        if [[ -z $extension ]]; then
            echo "no string extension specified for Sampling_Site"
            paste -d'\t' <(grep "analysis_project_id"  $config | cut -f2 -d' ') \
            <(grep "taxon_name" $config | rev | cut -f1 -d' ' | \
            rev | awk '{gsub(/^ +| +$/,"")} {print $0}') >> $gold_map
        else
            echo "string extension '${extension}' will be appended to Sampling_Site"
            paste -d'\t' <(grep "analysis_project_id"  $config | cut -f2 -d' ') \
            <(grep "taxon_name" $config | rev | cut -f1 -d' ' | \
            rev | awk '{gsub(/^ +| +$/,"")} {print $0}' | sed -e "s/$/$extension/") >> $gold_map
        fi

    else
        echo "WARNING: file $config does not exist, or has size 0. Ignoring."
    fi
done

if [[ $(wc -l <${gold_map}) -ge 2 ]]; then
    echo "Map successfully written to: $gold_map"
else
    echo "WARNING: no map records written to $gold_map"
fi
