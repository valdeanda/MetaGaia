# MetaGaia

*Sahil B. Shah - sbs2756@utexas.com, Valerie de Anda - valdeanda@utexas.edu*

# Project overview

MetaGaia will allow for the integrative analysis of both metagenomes and bins to better understand microbial communities.

The Gaia hypothesis, also known as Gaia theory or Gaia principle, proposes that all organisms and their inorganic surroundings on Earth are closely integrated to form a single and self-regulating complex system, maintaining the conditions for life on the planet

[Lovelock, J.E and Margulis, L. (1974). “Atmospheric homeostasis by and for the biosphere- The Gaia hypothesis”. Tellus 26 (1): 2–10](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.2153-3490.1974.tb01946.x)

## Requirements

This protocol assumes the following:
  
1. All necessary fastq files have been downloaded and placed within the same directory.
2. All necessary fna files have been downloaded and placed within the same directory.
3. A file has been created that maps each bin to its respective sample name (see "Mapping bins to samples" for more information).
4. The GTDBTK tool has been downloaded to map each bin to its repective taxa (see "Mapping taxa" for more information).
                                                                                                                    
# Usage notes

The goal of this program is to provide an end-to-end solution for calculating and visualizing the relative bin abundance by taking genome size and total read counts for each sample into account.

All output files will be saved in the `output` directory within the repository.

# Compute and visualize relative bin abundances

## 1. Creating bin abundance input files

1. To create all the input files necessary, run the `files_prep.py` script. Two input files are necessary and the directories where the `fastq` and `fna ` files are located are needed to run this script. These are the arguments necessary to run this script:

```
usage: files_prep.py [-h] -b BIN_SAMPLES -d DEPTH -q FASTQ_DIR -a FNA_DIR

optional arguments:
  -h, --help            show this help message and exit
  -b BIN_SAMPLES, --bin_samples BIN_SAMPLES
                        Input file path containing each bin mapped to each
                        sample with extension. Only tsv or csv files allowed
                        (str).
  -d DEPTH, --depth DEPTH
                        Input file path containing depth information (with
                        extension). Only tsv or csv files allowed (str).
  -q FASTQ_DIR, --fastq_dir FASTQ_DIR
                        Input directory path containing all the fasta files
                        (str).
  -a FNA_DIR, --fna_dir FNA_DIR
                        Input directory path containing all the fna files
                        (str).
```

To create the `bin_sample` file, read the "Mapping bins to samples" section below. To create the `depth` file, read the "Create initial depth file" section below as well.

Here is an example of how to run this script:

```
python3 files_prep.py -b bin2sample.tsv -d depth.tsv -q ~/fastq_files/ -a ~/fna_files/
```

2. There should be four files named "reads_file.tsv", "mapping_file.tsv", "depth_file", and "binsize_file.tsv" saved in the `output` folder.

## 2. Calculate relative bin abundance

1. The `bin_abundance.py` script then takes in a list of input files (processed from by running the `files_prep.py` script) that need to be created beforehand. All input files should be placed in the `data` directory. Here are the descriptions of each of the input files:

```
usage: bin_abundance.py [-h] -r READS -m MAPPING -d DEPTH -s BINSIZE
                        [-n READABLENUM]

optional arguments:
  -h, --help            show this help message and exit
  -r READS, --reads READS
                        Total number of reads. Rows are samples column are
                        reads
  -m MAPPING, --mapping MAPPING
                        Tabular file containingOriginal_Contig_Name Bin Sample
  -d DEPTH, --depth DEPTH
                        Tabular file depth infoOriginal_Contig_Name contigLen
                        Saple_Depth Depth
  -s BINSIZE, --binsize BINSIZE
                        Tabular file with Bin and corresponding Genome size
                        (bp) as columns
  -n READABLENUM, --readablenum READABLENUM
                        A large number to multiply the relative abundances so
                        that it is human readable.
```

Here is an example of how to run the script (using example data provided in the `data` folder):

```
python3 bin_abundance.py -r ../../data/example_input_files/example__reads.tsv -m ../../data/example_input_files/example__mapping.tsv -d ../../data/example_input_files/example__depth.tsv -s ../../data/example_input_files/example__size.tsv 1000000
```

2. The output will be a file saved in the `output` directory. This file can then be used in the next script: `bin_abundance_viz.py` to visualize the data.

## 3. Visualize relative bin abundance

1. With the `bin_abundance_viz.py` script, the files outputted from `bin_abundance.py` can be used to visualize the relative bin abundances in addition to a file mapping each bin to their respective taxa (see "Mapping taxa" section). Place all input files in the `data` directory. Here is an example of how to run the visualization script:

```
usage: bin_abundance_viz.py [-h] -b BIN_ABUNDANCE -t TAXONOMY_INFO
                            [-p PERCENT] [-w WIDTH] [-l HEIGHT] [-d DPI]
                            [-o OUT_FIG] [-c TAXA_COLOR]

optional arguments:
  -h, --help            show this help message and exit
  -b BIN_ABUNDANCE, --bin_abundance BIN_ABUNDANCE
                        Input file path outputted from MetaGaia with bin
                        abundances with extension (str).
  -t TAXONOMY_INFO, --taxonomy_info TAXONOMY_INFO
                        Input file path mapping taxonomy to each bin (with
                        extension).
  -p PERCENT, --percent PERCENT
                        Percent of highest sample in each bin [10] (float).
  -w WIDTH, --width WIDTH
                        Width of outputted clustermap figure [4] (int).
  -l HEIGHT, --height HEIGHT
                        Height of outputted clustermap figure [5] (int).
  -d DPI, --dpi DPI     Resolution for output figure file [300] (int).
  -o OUT_FIG, --out_fig OUT_FIG
                        Stores the figure in the specified file path and
                        format [test.png] (str).
  -c TAXA_COLOR, --taxa_color TAXA_COLOR
                        Input file path containing the RGB color code for each
                        taxa with extension [""] (str).
```

Here is an example of how to run this script:

```
python3 bin_abundance_viz.py -b ../output/bin_abundance_output.tsv -t ~/tax2bin.tsv -p 10 -w 8 -l 5 -d 500 -o figure1.jpg -c ~/tax2color.tsv
```

2. Outputs will be saved in the `output` directory containing the user-specified name.

# Additional commands possibly needed

## Mapping bins to samples

The user must manually map each bin to its respective sample.

## Create initial depth file

The script, `jgi_summarize_bam_contig_depths.pl`, must be used to create the depth file needed in `files_prep.py` and can be found here along with its usage notes: https://bitbucket.org/berkeleylab/metabat/src/master/.

## Mapping taxa

One way to map each the taxa to each bin is to use the GTDBTK tool found here: https://github.com/Ecogenomics/GTDBTk.