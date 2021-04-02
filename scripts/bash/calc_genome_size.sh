#!/bin/bash

cd $1

for i in *.fna; do grep -v ">" $i | tr -d '\n' | wc -c > $i.totalbp.txt; done
grep -vI “\x00” -- *.txt | sed 's/:/\t/g' > genomesize.tab
sed -i 's/.fna.totalbp.txt//g' genomesize.tab
