cd ../../data/IMG_data/

for i in *.tsv; do sed -i 's/ /_/g'  $i ;done 
#Edit (Sample_name) to actual sample name and repeat this line for each sample in dataset
sed -i '/s/scaffold/(Sample_name)_scaffold/g' the_mapping_file.tsv
for i in *.tsv ; do  awk -F "\t"  '{ print $3, $2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$21,$22,$23 }' $i > $i.sorted.tsv ; done
for i in *.sorted.tsv; do sed -i 's/ /\t/g' $i ; done
for i in *.sorted.tsv; do cut -f 1 $i  | sort | uniq   | sed 's/Original_Contig_Name//g' > $i.unique.scaffolds.tab ; done
cat *.unique.scaffolds.tab > all.scaffolds.tab