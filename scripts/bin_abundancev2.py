
#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     bin_abundancev2.py
# Purpose:  calculate the abundance of bins based of scaffold depth,
#           bin size (bp)
# @uthors:      vda - valdeanda@utexas.edu
#
# Created:     December 2020
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------


import argparse
import pandas as pd


# options
epilog = """Example:

$ python3 bin_abundancev2.py -r 20210104__reads.tab -m 20210104__mapping.tab -d 20210104__depth.tsv -s 20210104__size.tab"""

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

# Map the total reads of each sample  in the mapping file
merged1=pd.merge(df_mapping,df_reads[['Sample', 'Reads']], on='Sample')
# Map the Bin size to the mapping file
merged2=pd.merge(merged1,df_size[['Bin','Size']],on ='Bin')
# Map the coverage of each scaffold in different samples to the mapping file
merged3=pd.merge(merged2,df_depth[['Original_Contig_Name', 'contigLen', 'Sample_Depth', 'Depth']],on='Original_Contig_Name')

#Rename columns to perform a final merge

df=merged3.rename(columns={"Sample": "Sampling_Site","Sample_Depth": "Sample",'Reads':'Total_Reads'})
#Final dataframe to compute relative abundance
df=pd.merge(df,df_reads[['Sample', 'Reads']], on='Sample')
#Columns of the final data frame

#Original_Contig_Name	Bin	Sampling_Site	Total_Reads	Size	contigLen	Sample	Depth	Reads

# Calculate coverage
df['cov']=(df['contigLen'] *df['Depth'])/df['contigLen']

# Sum the coverage values by Bin
df['Sum_cov'] = df['cov'].groupby(merged3['Bin']).transform('sum')

# Normalize coverage values by Bin size
df['NormalizedCov'] = df['Sum_cov']/df['Size']

#Normalize by total number of reads of the mapped sample
df['RelativeAbundance'] = df['NormalizedCov']/df['Reads']

#Multiply by a big number to make the abundance readable
df['RelativeAbundanceReadable'] = df['RelativeAbundance']*100000000

#Drop duplicates

finaldf=df[['Bin', 'Sample',
       'NormalizedCov', 'RelativeAbundance', 'RelativeAbundanceReadable']]
finaldf=finaldf.drop_duplicates()

#Output files


mappingfile = args.mapping + "_intermediate_abundance_values.tsv"
abundance = args.mapping + "_abundanceby_bin.tsv"

df.to_csv(mappingfile,sep='\t')
finaldf.to_csv(abundance,sep='\t')


print("[END] Done computing abundance....................\n"
      "Please check the output files:\n",
      "Mapping File with intermediate values to compute abundance:", args.mapping + "_intermediate_abundance_values.tsv\n"
      "File containing relative abundance normalized by genome size and total read counts per sample:", args.mapping + "_abundanceby_bin.tsv\n"
      "Have a nice day :D")

