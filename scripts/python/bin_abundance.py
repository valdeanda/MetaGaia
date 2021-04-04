
#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     bin_abundance.py
# Purpose:  calculate the abundance of bins based of scaffold depth,
#           bin size (bp)
# @uthors:     vda  valdeanda@utexas.edu
# @ modif:     sbs sahil2699@gmail.com 
#              xgong xianzhe.gong@gmail.com
# Created:     December 2020
# Last updated Feb  2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------


import argparse
import pandas as pd

def main()

  # options
  epilog = """Example:

  $ python3 bin_abundancev3.py -r example__reads.tab -m example__mapping.tab -d example__depth.tsv -s example__size.tab"""

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

  args = parser.parse_args()


  reads        = args.reads
  mapping      = args.mapping
  depth        = args.depth
  size         = args.binsize

  df_reads    = pd.read_csv(reads,sep="\t")
  df_mapping  = pd.read_csv(mapping, sep="\t",index_col=False)
  df_size     = pd.read_csv(size, sep ="\t")
  df_depth    = pd.read_csv(depth, sep ="\t")


  #The user can change this number
  readablenum = 100000000

  # Map the total reads of each sample  in the mapping file

  #merged1=pd.merge(df_mapping,df_reads[['Sample', 'Reads']], on='Sample')
  # Map the Bin size to the mapping file
  mapping_size=pd.merge(df_mapping,df_size[['Bin','Size']],on ='Bin')
  # Map the coverage of each scaffold in different samples to the mapping file
  mapping_size_depth=pd.merge(mapping_size,df_depth[['Original_Contig_Name', 'contigLen', 'Sample_Depth', 'Depth']],on='Original_Contig_Name')

  #Rename columns to perform a final merge

  mapping_size_depth_new=mapping_size_depth.rename(columns={"Sample": "Sampling_Site","Sample_Depth": "Sample",'Reads':'Total_Reads'})
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

  finaldf=mapping_size_depth_new_reads[['Bin', 'Sample',
         'NormalizedCov', 'RelativeAbundance', 'RelativeAbundanceReadable']]
  finaldf=finaldf.drop_duplicates()

  #Output files

  mappingfile = args.mapping + "_IMGap_OUT_intermediate_abundance_values.tsv"
  abundance = args.mapping + "_IMGap_OUT_abundanceby_bin.tsv"

  mapping_size_depth_new_reads.to_csv("../../output/" + mappingfile, sep = '\t', index = False)
  finaldf.to_csv("../../output/" + abundance, sep = '\t', index = False)


  print("\n"
        "[END] Done computing abundance..........................................................................:\n"
        "Please check the output files:..........................................................................:\n\n"

        "1. Mapping File with intermediate values:", args.mapping + "__IMGap_OUT_intermediate_abundance_values.tsv\n"

        "2. File containing relative abundance normalized by bin:", args.mapping + "_IMGap_OUT_abundanceby_bin.tsv\n\n"

        "Have a nice day :D\n")

if __name__ == "__main__":
  main()

