#!/usr/bin/env bash

#Modified from Fasta_to_Scaffolds2Bin.sh:
#https://github.com/cmks/DAS_Tool/blob/master/src/Fasta_to_Scaffolds2Bin.sh

usage="
$(basename "$0"): create the contig to bin mapping file reqired for bin_abundance.py

where:

   -h --- show this help message
   -g --- directory containing genome FASTA
   -e --- extension of genome FASTA file(s) - must include full extension - e.g. .fna , NOT fna
    "

while getopts ':hg:e:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    g) genome_dir=${OPTARG};;
    #o) map_file=${OPTARG};;
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
test -d $genome_dir || { echo "ERROR: genome directory $genome_dir does not exist, exiting..."; exit 1; }

#printf "Original_Contig_Name\tBin\tSampling_Site\n" > $map_file

for genome in $genome_dir\/*${extension}; do
    bin_name=$(basename -s $extension $genome)
    grep ">" $genome | perl -pe "s/\n/\t$bin_name\n/g" | perl -pe "s/>//g"
done
