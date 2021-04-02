#!/bin/bash
  
fqfile=$1
nreads=$(echo "$(zcat $fqfile | wc -l)/4" | bc)
echo "$fqfile $nreads"
