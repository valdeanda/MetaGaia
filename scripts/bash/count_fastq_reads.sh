#!/bin/bash
for i in ../../fastq/*gz; do
  time zcat $i | perl ../JGI_tools/fastx-length.pl > /dev/null
  echo $i
done
