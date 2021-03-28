Bin_abun_prep.py

Use grep > to get all info from fna for mapping

User provides file mapping Bin to sample

Depth file derived from jgi_summarize_bam_contig_depths.pl script

----------------------------------------------------------------------

Add warning to count fastq script - Sample and Reads
Bin name matches bin name provided in bin-sample file (input)






# MetaGaia overview

MetaGaia will allow for the integrative analysis of both metagenomes and bins to better understand microbial communities. (vague)

The Gaia hypothesis, also known as Gaia theory or Gaia principle, proposes that all organisms and their inorganic surroundings on Earth are closely integrated to form a single and self-regulating complex system, maintaining the conditions for life on the planet

[Lovelock, J.E and Margulis, L. (1974). “Atmospheric homeostasis by and for the biosphere- The Gaia hypothesis”. Tellus 26 (1): 2–10](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.2153-3490.1974.tb01946.x)

## Requirements

This protocol assumes the following:
  
1. You have an account at the Integrated Microbial Genomes with Microbiome Samples system [IMG](https://img.jgi.doe.gov/cgi-bin/mer/main.cgi) 
2. You submitted your sequencing data  at the IMG system and waited for the annotation to be complete in all your samples
3. You have Metagenome Reconstructed Genomes (MAGs) and their corresponding genome size.
4. You will also need each MAG's taxonomic affiliation (this step can be done with automatic methods such as [GTDBk](https://github.com/Ecogenomics/GTDBTk), phylogeny using a set of conserved marker genes (i.e [Phylosift](https://github.com/gjospin/PhyloSift)), or 16S rRNA phylogenies. At the Baker Lab, we use a combination of the previously mentioned methods to assign taxonomy of MAGs. 
                                                                                                                    
# Usage notes

The goal of this program is to provide an end-to-end solution for calculating and visualizing the relative bin abundance by taking genome size and total read counts for each sample into account.

All output files will be saved in the `output` directory within the repository.

All bin abundance input files should be within the `data` directory.

# Compute and visualize relative bin abundances

## Calculate relative bin abundance

1. Obtain and prepare data files from IMG (see "Obtaining IMG data" section).

2. The `bin_abundance.py` script then takes in a list of input files (processed from IMG data) that need to be created beforehand. All input files should be placed in the `data` directory. Here are the descriptions of each of the input files:

```
python3 bin_abundance.py -h

usage: bin_abundance.py [-h] -r READS -m MAPPING -d DEPTH -s BINSIZE

optional arguments:
  -h, --help            show this help message and exit
  -r READS, --reads READS
                        Total number of reads. Rows are samples column are
                        reads
  -m MAPPING, --mapping MAPPING
                        Tabular file containing Original_Contig_Name Bin Sample
  -d DEPTH, --depth DEPTH
                        Tabular file depth info containing Original_Contig_Name Saple_Depth
                        Depth
  -s BINSIZE, --binsize BINSIZE
                        Tabular file with Bin and corresponding Genome size
                        (bp) as columns
```

Here is an example of how to run the script (using example data provided in the `data` folder):

```
python3 bin_abundance.py -r ../data/example_input_files/example__reads.tsv -m ../data/example_input_files/example__mapping.tsv -d ../data/example_input_files/example__depth.tsv -s ../data/example_input_files/example__size.tsv
```

2. Outputs will be saved in the `output` directory. Two files will be produced: an intermediate and a final (most important).

## Visualize relative bin abundance

1. With the `metagaia_viz.py` script, the files outputted from `bin_abundance.py` can be used to visualize the relative bin abundances in addition to a file mapping each bin to their respective taxa (see "Mapping taxa" section). Place all input files in the `data` directory. Here is an example of how to run the visualization script:

```
python3 [bin abundance output file with extension (string)] [bin taxonomy file with extension (string)] [top percent of samples within each bin (float)] [width of the outputted figures (int)] [height of the outputted figures (int)] [figure dpi (int)] [name of the outputted figures with extenstion (string)] [file with each taxa mapped to colors for visualization with extension (string)]
```

2. Outputs will be saved in the `output` directory containing the user-specified name.
                                                                                                                    
# Creating input files

To create all the input files necessary, run the `bin_abundance_prep.py` script.



Before creating the input files, complete the "Create file containing all the information" section first.

## Reads

To create the reads file, make sure all `.fastq files` are placed in the `fastq` directory and then run the `count_fastq_reads.sh` script like such:

```
sh count_fastq_reads.sh
```

Show example of file output.

Since the "Sample" names are paths, we need to edit them down to contain only the sample name. To do this, manually edit the name using the text editor of your choice (script is currently being developed).

## Mapping

To create the mapping file, run the `create_mapping.sh` script.

```
sh create_mapping.sh
```

Show example of file output.

## Depth

1. To create the depth file, run the `create_depth_file.sh` script.

```
sh create_depth_file.sh
```

Show example of file output.

2. Run the `modify_depth_file.py` script to format the file needed for calculating the relatvie bin abundance.

```
python3 modify_depth_file.py
```

Show example of file output.

## Binsize

To create the genome size file, make sure all `.fna files` are placed in the `fna` directory and then run the `calc_genome_size.sh` script like such:

```
sh calc_genome_size.sh
```

Show example of file output.

























# Create file containing all the information

This file contains all the columns needed for both the mapping and depth input files. This method consists of pulling the scaffolds from the `.tsv` consolidated file from IMG. Note, the names of the files must be changed to match the name of the sample for this method to work. To create this file, it involves a series of steps that results in:

Show example of file.

After changing the names of each of the IMG files, the steps needed to create this file are as follows:

1. Edit the `mapping_scaffolds.sh` script by adding the sample in the scaffold name. Thus, this needs to be done individually by copying and pasting this line with differing sample names:

```
sed -i '/s/scaffold/(Sample_name)_scaffold/g' the_mapping_file.tsv
```

2. Run the `mapping_scaffolds.sh` script to create a list of all the scaffolds of all your samples. The general name format within the list is `sample_name_scaffold_name`.

3. Create a file with the scaffolds that are binned by running the `xyz.sh` script.

4. Run the `mapping_scaffolds2.sh` script to concatenate the two scaffold files and map scaffolds without a bin.

Example:

```
sh mapping_scaffolds2.sh all_scaffolds.tab scaffold2bins.tab
```

5. Run the `get_len_gc.sh` script to obtain the length and GC content from each sample assembly file `.fa`. `scaffold2gclength.tab` will be outputted.

6. Run the `mapping_awk.sh` script to concatenate the `mappingFile1.tab` and `scaffold2gclength.tab ` files. The command should look like this:

```
sh mapping_awk.sh mappingFile1.tab scaffold2gclength.tab
```

7. Need help with getting depth info.

To be continued...

# Additional commands possibly needed

## Obtaining IMG data

Need help with this section since I have not gone through the process myself.

## Mapping taxa

Should I include these steps?
















