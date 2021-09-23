#!/usr/bin/env bash

#Author: Ian Rambo
#Created: September 22, 2021
#Purpose: Filter VirSorter2 and CheckV hits and extract passing viral genomes.

checkv_summary=$1
virsorter_summary=$2
virsorter_fna=$3
output_dir=$4

has_command () {
    #Make sure you can execute certain commands
    command -v "$1" >/dev/null 2>&1 || { echo "Requires $1. Ensure that $1 is in your \$PATH."; exit 1; }
}

has_command pullseq

usage="
$(basename "$0"): Quality screen VirSorter2 and CheckV hits and extract viral genomes.

where:

   -h --- show this help message
   -c --- path for CheckV quality_summary.tsv file
   -v --- path for VirSorter2 final-viral-score.tsv file
   -f --- path for VirSorter2 final-viral-combined.fa FASTA file
   -o --- output directory

   Screening metrics are described here: https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-bwm5pc86
    "

while getopts ':hi:o:e:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    c) checkv_summary=${OPTARG};;
    v) virsorter_summary=${OPTARG};;
    f) virsorter_fna=${OPTARG};;
    o) output_dir=${OPTARG};;

    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))

test -d $output_dir || mkdir -p $output_dir

keep_1=${output_dir}/uvig_keep_1.fna
keep_2=${output_dir}/uvig_keep_2.fna
combined=${output_dir}/uvig_combined.fna
# keep_2a=${output_dir}/uvig_keep_2a.fna
# keep_2b=${output_dir}/uvig_keep_2b.fna



echo "extracting UViGs keep_1: one or more viral gene, high viral probability"
#viral gene > 0 - CheckV should call these viral.
tail -n +2 $checkv_summary | awk '$6 > 0' | \
    pullseq -N -i $virsorter_fna > $keep_1

echo "wrote UViGs keep_1"


#keep_2
cat <( tail -n +2 $checkv_summary | awk '{ if( $6==0 && $7==0) print $0 }' | cut -f1 | sort -u ) \
    <( grep -f <(tail -n +2 $checkv_summary | awk '{ if( $6==0 ) print $0 }' | \
         cut -f1 | sort -u) <(tail -n +2 $virsorter_summary) | \
         awk '{ if( $5>=0.95 || $8>2 ) print $1 }' | \
         sort -u ) | \
         sort -u | \
         pullseq -N -i $virsorter_fna > $keep_2

cat $keep1 $keep2 > $combined
# echo "extracting UViGs keep_2a: no viral gene, no host gene"
# #No viral gene, no host gene
# tail -n +2 $checkv_summary | awk '{ if( $6==0 && $7==0) print $0 }' | \
#      cut -f1 | \
#      sort -u | \
#      pullseq -N -i $virsorter_fna > $keep_2a
# echo "wrote UViGs keep 2a"
#
# echo "extracting UViGs keep_2b: viral gene, high VirSorter2 score, 2 or more hallmark genes"
# #No viral gene, high VirSorter2 score, 2 or more hallmark genes
#  grep -f <(tail -n +2 $checkv_summary | awk '{ if( $6==0 ) print $0 }' | \
#      cut -f1 | sort -u) <(tail -n +2 $virsorter_summary) | \
#      awk '{ if( $5>=0.95 || $8>2 ) print $1 }' | \
#      sort -u | \
#      pullseq -N -i $virsorter_fna > $keep_2b && \
#      echo "wrote UViGs keep_2b"

echo "DONE"
