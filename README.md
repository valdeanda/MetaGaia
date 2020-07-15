# IMG_annotation protocol




---

#### Required files for bin_abundace.py  script

1. **Total lenght of the metagenomic assemblies in bp**

Generate a file with the total lenght of your assemblies (i.e lenghtpersample.tab) 

```bash
for i in *.fasta ; do grep -v ">" $i | tr -d '\n' | wc -c > $i.totalbp.txt ; done

#grep -v ">" find all non header lines (i.e. without ">")
#tr -d '\n' delete all newlines
#wc -c count all characters
```
Show the content of all your files and redirect the output to a new one

```bash
grep -vI "\x00" -- *.txt | sed 's/:/\t/g' > lenghtpersample.tab
```

This is an example of how the final should look like. Keep in mind that the name of the columns have to match the example files.

```
Sample  total_length(bp)
4484    55507215
4572    516163358
AB_03   896120870
AB_1215 700752333
```
2. **Tabular file containing original contig name Bin, GC, length and Depth in columns**

Extensive details of how to create this file soon

```
ID	Original_Contig_Name	Bin	GC	Length	Depth	Sample
0	AB_1215_scaffold_0	NoBin	0.315	219188	27.6207	AB_1215
1	AB_1215_scaffold_1	NoBin	0.348	203294	36.76	AB_1215
...
8	AB_1215_scaffold_100002	AB_1215_Bin_149	0.389	2338	14.8501	AB_1215
```
3.  **Tabular file contaning manually curated taxonomy according to phylogeny**

Extensive details of how to create this file soon. Complete file provided in Supplementary Material De Anda et al., In prep 

```
Bin     Domain  Class   Proposed taxonomy       Superphylum/Phylum      Notes   
AB_03_Bin_74    Bacteria        Acidobacteria   Acidobacteria bacterium AB_03_Bin_74    Acidobacteria   Taxonomy Confirmed with Phylogeny
AB_03_Bin_78    Bacteria        Acidobacteria   Acidobacteria bacterium AB_03_Bin_78    Acidobacteria   Taxonomy Confirmed with Phylogeny
AB_03_Bin_79    Bacteria        Acidobacteria   Acidobacteria bacterium AB_03_Bin_79    Acidobacteria   Taxonomy Confirmed with Phylogeny
AB_1215_Bin_144 Bacteria        Acidobacteria   Acidobacteria bacterium AB_1215_Bin_144 Acidobacteria   Taxonomy Confirmed with Phylogeny
```

# Analyzing the outputs 

## Bin abundance

```bash
python3 bin_abundance.py -l lenghtpersample.tab -m scaffolds2sample.tab -t Taxonomy.tsv
```
Get the taxonomy of the bins from the output file 

```bash
 cut -f 4,17 Taxonomy.tsv_abundance.tsv   | sort | uniq | grep -v "NoBin" > Taxonomy.tsv_abundance_bins.tsv 
 ```
