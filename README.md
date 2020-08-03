# IMGap
IMG **a**nnotation **p**rotocol Baker Lab

This protocol assumes the following:

1. You have an account at the Integrated Microbial Genomes with Microbiome Samples system [IMG](https://img.jgi.doe.gov/cgi-bin/mer/main.cgi) 
2. You submitted your sequencing data  at the IMG system and waited for the annotation to be complete in all your samples
3. You have Metagenome Reconstructed Genomes (MAGs) and their corresponding taxonomic affiliations (this step can be done with automatic methods such as [GTDBk](https://github.com/Ecogenomics/GTDBTk) or phylogenies using a set of conserved marker genes (i.e [Phylosift](https://github.com/gjospin/PhyloSift)) or 16S rRNA phylogenies. At the Baker Lab, we use a combination of the previosly mentioned methods to assign taxonomy of MAGs. 



## Download

Annotate metagenomes at the [IMG JGI platafform](https://img.jgi.doe.gov/cgi-bin/mer/main.cgi) 

After completion, you will download a directory of each of your metagenomic samples containing the following subdirectories: 

1. **IMG Data**: Contanining the following files 

```
207502.assembled.COG
207502.assembled.crisprs
207502.assembled.EC
207502.assembled.faa
207502.assembled.fna
207502.assembled.gff
207502.assembled.KO
207502.assembled.last.blasttab
207502.assembled.names_map
207502.assembled.Pfam.blout
207502.assembled.phylodist
207502.assembled.product_names
207502.assembled.unsorted.Pfam.hmmout
3300033388
```

2. **Metagenome Report Tables**: containing the following files 

```
Table_6_-_207017.annotation_parameters.txt
Table_7_-_3300033141.metagenome_properties.txt
Table_8_-_3300033141.taxonomic_composition.txt
Table_9_-_3300033141.functional_diversity.txt
```

## Consolidate data 

Using the script consolidateJGIdata.pl previously available at [JGI tools scripts](https://github.com/Geo-omics/scripts/tree/master/JGITools). This script is required  to create a single annotation table derived from JGI 
```{bash, highlight=TRUE, eval=FALSE}

 perl scripts/consolidateJGIdata.pl -DIR 3300033141/ -OUTDIR test -scripts /home/Guaymas_C_project/scripts/JGI_tools/
# consolidateJGIdata.pl v0.1.6
[LGC] Calculating length and GC content of all contigs...
[COG] File found:
[PFAM] File found:	3300033141/3300033141.a.pfam.txt
[TAXA] File found:	3300033141/3300033141.a.phylodist.txt
[KO] File found:	3300033141/3300033141.a.ko.txt
[EC] File found:	3300033141/3300033141.a.ec.txt
[PROD] File found:	3300033141/3300033141.a.gene_product.txt
[MAP] File found:	3300033141/3300033141.a.map.txt

# Bin	NumContigs	Fraction_of_Dataset	NumGenes	Fraction_of_Dataset	Total_Length_of_Bin	Average_GC	Number_of_CDS
Unclassified	211043	1	1419882	1	0	0	1400309

```

If everything goes right,  a tsv  file will be created.  That file contains the mapping scaffold names and annotation










---

# Required Files

## Bin abundance in metagenomic samples

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

Extensive details of how to create this file below, section: Creating mapping Files

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


## Bin abundance

```bash
python3 bin_abundance.py -l lenghtpersample.tab -m scaffolds2sample.tab -t Taxonomy.tsv
```
Get the taxonomy of the bins from the output file 

```bash
 cut -f 3,16 Taxonomy.tsv_abundance.tsv   | sort | uniq | grep -v "NoBin" > Taxonomy.tsv_abundance_bins.tsv 
 ```
 
 
# Creating Mapping Files 

1. You need a file containing the list of all the scaffolds of all your samples. For purpuses of this example this file will be named *all_scaffolds.tab*, you can give it the name that you want. You can create this file 

```bash
less all_scaffolds.tab

sample0_scaffold_0
sample1_scaffold_1
sample2_scaffold_10
sample3_scaffold_100
sample4_scaffold_1000
```
To create this file, there are several options, one way to do it, is to pull of the scaffolds from the TSV consolidate file from IMG but you need to change the name first to make sure that  matches the name of the sample at the beggining of each scaffold name. Here is one way to do it.  

1.1. Delete white spaces:
```bash
for i in *.tsv; do sed -i 's/ /_/g'  $i ;done 
```


1.2. Add the name of the sample in the scaffold name. This need to be done for all samples individually.

```bash
sed -i '/s/scaffold/(Sample_name)_scaffold/g' the_mapping_file.tsv

```


1.3. Sort file so the first column are the scaffolds  and delete white columns 

```bash
for i in *.tsv ; do  awk -F "\t"  '{ print $3, $2,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$21,$22,$23 }' $i > $i.sorted.tsv ; done
```


1.4. Substitute spaces by tab  

```{bash, highlight=TRUE, eval=FALSE}
for i in *.sorted.tsv; do sed -i 's/ /\t/g' $i ; done
```


1.5. Get the scaffolds for the mapping file. 

```bash
for i in *.sorted.tsv; do cut -f 1 $i  | sort | uniq   | sed 's/Original_Contig_Name//g' > $i.unique.scaffolds.tab ; done 

cat *.unique.scaffolds.tab > all.scaffolds.tab
```


2. Then, you need to create a file  file with the scaffold that are binned, in this example *scaffolds2bins.tab*
Note that scaffold name and the bin have to contain the name of the sample of the scaffold name

```bash
less scaffolds2bins.tab
sample1_scaffold_1 sample1_Bin_7
sample2_scaffold_10 sample2_Bin_14
sample3_scaffold_10003 sample3_Bin_27
sample4_scaffold_10005 sample4_Bin_9
```

3. Create a file that contains the scaffolds that are not binned. Combine  *all_scaffolds.tab* and *scaffolds2bins.tab* 

```bash

bash scripts/mapping_awk.sh all_scaffolds.tab scaffolds2bins.tab
#After using the mapping_awk.sh script remove ALWAYS remove matching characters and change spaces by tabs 
sed  -e "s/\r//g"  all_scaffolds.tab.mapping.tsv |  sed 's/ /\t/g' >  all_scaffolds.tab.mapping_sorted.tsv 
```

3.1 The above awk script will  create concatenate both files file all_scaffolds.tab.mapping_sorted.tsv 

```bash
less all_scaffolds.tab.mapping_sorted.tsv 

sample1_scaffold_0
sample1_scaffold_1
sample1_scaffold_1 sample1_Bin_7
sample2_scaffold_10
sample2_scaffold_10 sample2_Bin_14
sample2_scaffold_100
WB1_scaffold_1000
```
3.2 Now we need to specify which ones are not binned, in this case are the ones with empty spaces, and we are going to create *mappingFile1.tab*

```bash
cut -f 1,2 all_scaffolds.tab.mapping_sorted.tsv  | awk '{if (!$2) {print $1, "NoBin"} else {print $1, $2}}' > mappingFile1.tab 
#again, remove matching characters and replace by tabs

sed -i -e "s/\r//g" mappingFile1.tab 
sed -i 's/ /\t/g' mappingFile1.tab 

less mappingFile1.tab

sample1_scaffold_0 NoBin
sample1_scaffold_1 NoBin
sample1_scaffold_1 sample1_Bin_7
sample2_scaffold_10 NoBin
sample2_scaffold_10 sample2_Bin_14
sample2_scaffold_100 NoBin
WB1_scaffold_1000 NoBin
```

3. Now we need to compute the %GC and lenght(bp) of each scaffold and create a mapping file *scaffold2gclength.tab*

```bash
less scaffold2_gc_length.tab

sample1_scaffold_0  0.346   344315
sample1_scaffold_1  0.652   402571
sample1_scaffold_2  0.624   298614
sample2_scaffold_3  0.332   365716
sample2_scaffold_4  0.649   305259
sample2_scaffold_5  0.642   212965

```
3.1 This file is created with the [length+GC.pl](https://github.com/valdeanda/IMG_annotation/blob/master/JGI_tools/length%2BGC.pl) script provided in the JGI_tools directory 

The input is the metagenomic data that was submitted to JGI, and you have to do it for all the samples with a simple bash loop 

```bash
#For one sample
perl JGI_tools/length+GC.pl sample_assemlie.fa > scaffold2gclength.tab

#For all your assemblies located in a specific directory 

for sample in directory/*.fa do perl JGI_tools/length+GC.pl $sample > $sample.scaffold2length.tab ; done
cat *.scaffold2length.tab >  scaffold2gclength.tab 

```

4. Create the final mapping file *mappingFile1.tab.mapping.tsv*  that will contain scaffold information of bins, %GC and lenght. Use again the awk script to concatenate both files *mappingFile1.tab* and *scaffold2gclength.tab* 

```bash
./mapping_awk.sh mappingFile1.tab and scaffold2gclength.tab 

#The above command will create the mappingFile1.tab.mapping.tsv 

sample1_scaffold_0  NoBin   0.346   344315
sample1_scaffold_1  NoBin   sample1_Bin_7       0.652   402571
sample1_scaffold_10 NoBin   sample1_Bin_14      0.328   189441
sample1_scaffold_100        NoBin   0.384   96095

#again, remove matching characters and replace by tabs

sed -i -e "s/\r//g" mappingFile1.tab.mapping.tsv 
sed -i 's/ /\t/g' mappingFile1.tab.mapping.tsv  

```
5. The next File contains the Depth information of each scaffold

**More details of how to create the Depth information soon**

```
Original_Contig_Name	Depth
sample1_scaffold_0	20.3781
sample1_scaffold_100000	10.585
sammple1_scaffold_100001	49.4958
sample1_scaffold_100002	11.9636
```
5.1 Use the script joinFile.sh to add the depth information to the mappingFile

```bash
./joinFile.sh MappingFile4col.tab GBdepth.tab  > MappingFile5col.tab
```

5. Create the last mapping file that contains the name of your sample and the following 7 column headers. 
Please make sure that you have the exact name and no space after the name of the column.! 


```bash
sed 's/_scaffold/\t/g' MappingFile5col.tab | cut -f 1 | sed 's/Original_Contig_Name/Sample/g' > MappingSamples.txt
 paste MappingFile5col.tab MappingSamples.txt > MappingFile6col.tab

ID	Original_Contig_Name	Bin	GC	Length	Depth	Sample
0	AB_1215_scaffold_0	NoBin	0.315	219188	27.6207	AB_1215
1	AB_1215_scaffold_1	NoBin	0.348	203294	36.76	AB_1215
2	AB_1215_scaffold_10	NoBin	0.379	119163	50.7298	AB_1215
3	AB_1215_scaffold_100	NoBin	0.4	54468	31.3466	AB_1215
```


# Useful commands 

### Download fastq data from IMG 


Download from IMG <https://genome.jgi.doe.gov/portal/Sample.info.html>.

Follow instructions for IMG downloads, "Download with API": <https://genome.jgi.doe.gov/portal/help/download.jsf#/api>.

1. Log in to IMG account on server (enter password):

```{bash, eval = FALSE}
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=x@email.com' --data-urlencode 'password=x' -c cookies > /dev/null
```

2. Download a list of files available for the portal that you are interested in. Fill in your unique organism after "organism="

An example portal:
```{bash, eval = FALSE}
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism=X' -b cookies > files.xml
```

3. Find the fastq in the XML doc under "Filtered Raw Data". Paste the provided url into the command below, starting at "/portal/" for file "XXX-filter-METAGENOME.fastq.gz". ** In my experience the url provided on the XML doc online works, while the url in the XML doc downloaded to the server does not provide the correct link.

```{bash, eval = FALSE}
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/X.CCAAGCA-TTGCTTG.filter-METAGENOME.fastq.gz' -b cookies > X-METAGENOME.fastq.gz
```

Use Sickle to trim fastqs: If you have one file with interleaved reads as input and you want ONLY one interleaved file as output:

```{bash, eval = FALSE}
sickle pe -g -c X-METAGENOME.fastq.gz -t sanger -M X.fastq.gz
```
----

### Download assemblies data from IMG 


Replace X with your actual data 

1. Log onto IMG on server

```bash
curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=youraccount@gmail.com' --data-urlencode 'password=yourpassword!' -c cookies > /dev/null
```
2. Download tar.gz (link in XML file at <https://genome.jgi.doe.gov/portal/X.info.html>)

```bash
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/X/XXXXXXXXXX.tar.gz' -b cookies > XXXXXXXXXX.tar.gz
```
3. Uncompress file

```bash
tar xvzf XXXXXXXXXX.tar.gz
```

4. Download Assemblies

```bash 
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/X.contigs.fasta' -b cookies > XXXXX.contigs.fasta
```


