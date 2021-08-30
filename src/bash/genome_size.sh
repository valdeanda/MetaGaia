#!/bin/bash
#Author: Ian Rambo
#Created: August 25, 2021
#Purpose: Calculate genome sizes (bp) for a directory of genome FASTA files

usage="
$(basename "$0"): create the binsize file reqired for bin_abundance.py

where:

   -h --- show this help message
   -g --- directory containing genome FASTA
   -o --- output binsize file
   -e --- extension of genome FASTA file(s) - must include full extension - e.g. .fna , NOT fna
    "

while getopts ':hg:o:e:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    g) genome_dir=${OPTARG};;
    o) size_file=${OPTARG};;
    e) ext=${OPTARG};;

    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))
#-----------------------------------------------------------------------------
test -d $genome_dir || { echo "ERROR: genome directory $genome_dir does not exist, exiting..."; exit 1; }

out_root=$(dirname $size_file)
test -d $out_root || mkdir -p $out_root

size_file_tmp="${size_file}.tmp"

for genome in $genome_dir/*${ext}; do
    genome_size=$(grep -v '>' $genome | tr -d '\n' | wc -c)
    genome_name=$(basename -s $ext $genome)
    printf "$genome_name\t$genome_size\n"; done >${size_file_tmp}

printf "Bin\tSize\n" | cat - $size_file_tmp > ${size_file} && rm $size_file_tmp


echo "FINISHED!"
