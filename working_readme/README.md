Bin_abun_prep.py

Depth file derived from jgi_summarize_bam_contig_depths.pl script






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

## Creating bin abundance input files

To create all the input files necessary, run the `bin_abundance_prep.py` script. Two input files are necessary and the directories where the `fastq` and `fna ` files are needed to run this script. These are the arguments necessary to run this script:

```
usage: bin_abundance_prep.py [-h] bin_samples depth fastq_dir fna_dir

positional arguments:
  bin_samples  Input file path containing each bin mapped to each sample with
               extension. Only tsv or csv files allowed. (str)
  depth        Input file path containing depth information (with extension).
               Only tsv or csv files allowed. (str)
  fastq_dir    Input directory path containing all the fasta files. (str)
  fna_dir      Input directory path containing all the fna files. (str)
```

To create the `bin_sample` file, read the "Mapping taxa" section below. To create the `depth` file, read the "Format depth file" section below as well.

Here is an example of how to run this script:

```
python3 bin_abundance_prep.py bin2sample.tsv depth.tsv ~/fastq_files/ ~/fna_files/
```

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
python3 bin_abundance.py -r ../../data/example_input_files/example__reads.tsv -m ../../data/example_input_files/example__mapping.tsv -d ../../data/example_input_files/example__depth.tsv -s ../../data/example_input_files/example__size.tsv
```

2. Outputs will be saved in the `output` directory. Two files will be produced: an intermediate and a final (most important).

## Visualize relative bin abundance

1. With the `metagaia_viz.py` script, the files outputted from `bin_abundance.py` can be used to visualize the relative bin abundances in addition to a file mapping each bin to their respective taxa (see "Mapping taxa" section). Place all input files in the `data` directory. Here is an example of how to run the visualization script:

```
usage: metagaia_viz.py [-h]
                       bin_abundance taxonomy_info [percent] [width] [height]
                       [dpi] [out_fig] [taxa_color]

positional arguments:
  bin_abundance  Input file path outputted from MetaGaia with bin abundances
                 with extension (str).
  taxonomy_info  Input file path containing taxonomy information (with
                 extension).
  percent        Percent of highest sample in each bin [10] (float).
  width          Width of outputted clustermap figure [4] (int).
  height         Height of outputted clustermap figure [5] (int).
  dpi            Resolution for output figure file [300] (int).
  out_fig        Stores the figure in the specified file path and format
                 [test.png] (str).
  taxa_color     Input file path containing the color code for each taxa with
                 extension [""] (str).
```

Here is an example of how to run this script:

```
python3 metagaia_viz.py bin_abundance_output.tsv tax2bin.tsv 10 8 5 500 figure1.jpg tax2color.tsv
```

2. Outputs will be saved in the `output` directory containing the user-specified name.

# Additional commands possibly needed

## Obtaining IMG data

Need help with this section since I have not gone through the process myself.

## Mapping taxa

Should I include these steps?

## Format depth file
















