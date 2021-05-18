#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     bin_abundance.py
# Purpose:  calculate the abundance of bins based of scaffold depth,
#           bin size (bp)
# @uthors:     vda - valdeanda@utexas.edu
# @ modif:     sbs - sahil2699@gmail.com 
#              xgong xianzhe.gong@gmail.com
# Created:     December 2020
# Last updated Feb  2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------


import argparse
import os
import pandas as pd


def main():

  # options
  epilog = """Example:

  $ python3 bin_abundancev3.py -r example__reads.tsv -m example__mapping.tsv -d example__depth.tsv -s example__size.tsv"""

  parser = argparse.ArgumentParser(description=__doc__, epilog=epilog)
  parser.add_argument('-r','--reads', required=True,
                      help=("Total number of reads. Rows are samples "
                      "column are reads "))

  parser.add_argument('-m', '--mapping', required=True,
                      help=("Tabular file containing"
                            "Original_Contig_Name	Bin  Sample"))

  parser.add_argument('-d', '--depth', required=True,
                      help=("Tabular file depth info"
                          "Original_Contig_Name  contigLen  Saple_Depth Depth"))

  parser.add_argument('-s', '--binsize', required=True,
                      help=("Tabular file with Bin and corresponding Genome size (bp) as columns"))

  parser.add_argument('-n', '--readablenum', required=False, default=100000000,
                      help=("A large number to multiply the relative abundances so that it is human readable."))

  args = parser.parse_args()

  #Create output directory if not already present
  if not os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/../../output"):
    os.makedirs(os.path.dirname(os.path.abspath(__file__)) + "/../../output")


  reads        = args.reads
  mapping      = args.mapping
  depth        = args.depth
  size         = args.binsize
  readablenum  = float(args.readablenum)

  df_reads    = pd.read_csv(reads,sep="\t")
  df_mapping  = pd.read_csv(mapping, sep="\t",index_col=False)
  df_size     = pd.read_csv(size, sep ="\t")
  df_depth    = pd.read_csv(depth, sep ="\t")


  # Map the total reads of each sample  in the mapping file

  # Map the Bin size to the mapping file
  mapping_size=pd.merge(df_mapping,df_size[['Bin','Size']],on ='Bin')
  # Map the coverage of each scaffold in different samples to the mapping file
  mapping_size_depth_new=pd.merge(mapping_size,df_depth[['Original_Contig_Name', 'contigLen', 'Sample', 'Depth']],on='Original_Contig_Name')

  #Final dataframe to compute relative abundance

  mapping_size_depth_new_reads=pd.merge(mapping_size_depth_new,df_reads[['Sample', 'Reads']], on='Sample')
  #Columns of the final data frame

  #Original_Contig_Name	Bin	Sampling_Site	Total_Reads	Size	contigLen	Sample	Depth	Reads

  # Calculate coverage


  mapping_size_depth_new_reads['cov']=mapping_size_depth_new_reads['contigLen'].astype(float) * mapping_size_depth_new_reads['Depth'].astype(float)

  # Sum the coverage values by Bin

  mapping_size_depth_new_reads['Sum_cov'] = mapping_size_depth_new_reads.groupby(by=['Bin',"Sample"])['cov'].transform('sum')

  # Normalize coverage values by Bin size
  mapping_size_depth_new_reads['NormalizedCov'] = mapping_size_depth_new_reads['Sum_cov'].astype(float) / mapping_size_depth_new_reads['Size'].astype(float)

  #Normalize by total number of reads of the mapped sample
  mapping_size_depth_new_reads['RelativeAbundance'] = mapping_size_depth_new_reads['NormalizedCov'].astype(float) / mapping_size_depth_new_reads['Reads'].astype(float)

  #Multiply by a big number to make the abundance readable
  mapping_size_depth_new_reads['RelativeAbundanceReadable'] = mapping_size_depth_new_reads['RelativeAbundance'].astype(float) * readablenum

  #Drop duplicates
  finaldf=mapping_size_depth_new_reads[['Bin', 'Sample', 'Sampling_Site', 'RelativeAbundance', 'RelativeAbundanceReadable']]
  finaldf=finaldf.drop_duplicates()

  #Output files
  abundance = "MetaGaia_OUT_abundanceby_bin.tsv"

  finaldf.to_csv(os.path.dirname(os.path.abspath(__file__)) + "/../../output/" + abundance, sep = '\t', index = False)


  print("\n"
        "[END] Done computing abundance..........................................................................:\n"
        "Please check the output file:..........................................................................:\n\n"

        "1. File containing relative abundance normalized by bin: MetaGaia_OUT_abundanceby_bin.tsv\n\n"

        "Have a nice day :D\n")

if __name__ == "__main__":
  main()

