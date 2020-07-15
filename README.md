# IMG_annotation protocol Baker Lab

This protocol assumes the following:

1. You have an account at the Integrated Microbial Genomes with Microbiome Samples system [IMG](https://img.jgi.doe.gov/cgi-bin/mer/main.cgi) 
2. You submitted your sequencing data  at the IMG system and waited for the annotation to be complete in all your samples
3. You have Metagenome Reconstructed Genomes (MAGs) and their corresponding taxonomic affiliations (This step can be done with automatic methods such as [checkM](https://ecogenomics.github.io/CheckM/)or phylogenies using a set of conserved marker genes (i.e [Phylosift](https://github.com/gjospin/PhyloSift)) or 16S rRNA phylogenies. We use a combination of the previosly mentioned methods to assign taxonomy of MAGs. 



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
 
 
 
 
 

# Useful commands 

### Download fastq data from IMG 


Download from IMG <https://genome.jgi.doe.gov/portal/SpaChaFrBayDelta/SpaChaFrBayDelta.info.html>.

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


