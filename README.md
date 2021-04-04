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

# Compute and visualize relative bin abundances

## 1. Creating bin abundance input files

To create all the input files necessary, run the `files_prep.py` script. Two input files are necessary and the directories where the `fastq` and `fna ` files are located are needed to run this script. These are the arguments necessary to run this script:

```
python3 files_prep.py [a tsv or csv file mapping bins to their respective sample names (str)] [a tsv or csv file containing the depth information (str)] [the directory containing all the fastq files (str)] [the directory containing all the fna files (str)]
```

To create the `bin_sample` file, read the "Mapping bins to samples" section below. To create the `depth` file, read the "Create initial depth file" section below as well.

Here is an example of how to run this script:

```
python3 bin_abundance_prep.py bin2sample.tsv depth.tsv ~/fastq_files/ ~/fna_files/
```

## 2. Calculate relative bin abundance

1. Obtain and prepare data files from IMG (see "Obtaining IMG data" section).

2. The `bin_abundance.py` script then takes in a list of input files (processed from by running the `files_prep.py` script) that need to be created beforehand. All input files should be placed in the `data` directory. Here are the descriptions of each of the input files:

```
python3 bin_abundance.py [a tsv file containing the total number of reads in each sample (str)] [a tsv file mapping the contig name to its respective bin and sample names (str)] [a tsv file containing the depth information in a long tabular format (str)] [a tsv file mapping each bin with its corresponding genome size (str)]
```

Here is an example of how to run the script (using example data provided in the `data` folder):

```
python3 bin_abundance.py -r ../../data/example_input_files/example__reads.tsv -m ../../data/example_input_files/example__mapping.tsv -d ../../data/example_input_files/example__depth.tsv -s ../../data/example_input_files/example__size.tsv
```

2. Outputs will be saved in the `output` directory. Two files will be produced: an intermediate and a final (most important).

## 3. Visualize relative bin abundance

1. With the `metagaia_viz.py` script, the files outputted from `bin_abundance.py` can be used to visualize the relative bin abundances in addition to a file mapping each bin to their respective taxa (see "Mapping taxa" section). Place all input files in the `data` directory. Here is an example of how to run the visualization script:

```
python3 bin_abundance_viz.py [the final output file from the bin_abundance.py script (str)] [a tsv or csv file mapping each bin with its respective taxa (str)] [percent of the highest sample in each bin [10] (float)] [width of outputted clustermap [4] (int)] [height of outputted clustermap [5] (int)] [resolution of output figure [300] (int)] [desired name of the outputted figure with the extension [test.png] (str)] [a tsv or csv file containing the color code for each taxa [""] (str)]
```

Here is an example of how to run this script:

```
python3 metagaia_viz.py bin_abundance_output.tsv tax2bin.tsv 10 8 5 500 figure1.jpg tax2color.tsv
```

2. Outputs will be saved in the `output` directory containing the user-specified name.

# Additional commands possibly needed

## Mapping bins to samples

The user must manually map each bin to its respective sample.

## Create initial depth file

The script, `jgi_summarize_bam_contig_depths.pl`, must be used to create the depth file needed in `files_prep.py` and can be found here along with its usage notes: https://bitbucket.org/berkeleylab/metabat/src/master/.

## Mapping taxa

One way to map each the taxa to each bin is to use the GTDBTK tool found here: https://github.com/Ecogenomics/GTDBTk.