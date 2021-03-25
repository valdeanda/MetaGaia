for sample in ../../data/IMG_data/fa/*.fa do perl ../JGI_tools/length+GC.pl $sample > $sample.scaffold2length.tab ; done
cat *.scaffold2length.tab >  ../scaffold2gclength.tab 