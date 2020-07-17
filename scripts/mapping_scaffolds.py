#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     mapping_scaffolds.py
# Purpose:  use the IMG derived annotation file and map to the original scaffold names
#
# @uthors:      vda - valdeanda@utexas.edu, acph  dragopoot@gmail.com
#
# Created:     2019
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------

# import os
import argparse
# import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd


# options
epilog = """Example:

$ python3 mapping_scaffolds.py  -i 3300033141.tsv  -s mappingFile2.tab"""

parser = argparse.ArgumentParser(description=__doc__, epilog=epilog)
parser.add_argument('-i', '--filename', required=True,
                    help=("Input file in tabular format. Rows are genes "
                    "columns are IMG annotations."))

parser.add_argument('-s', '--scaffolds', required=True,
                    help=('Tabular file containin sample, bin'
                    'and original contig name columns'))

args = parser.parse_args()


fname    = args.filename
scaffold = args.scaffolds

df          = pd.read_csv(fname,index_col=0,sep="\t")
df_scaff    = pd.read_csv(scaffold,sep="\t")




# Map scaffolds 

df2=pd.merge(df,df_scaff[['Original_Contig_Name', 'Bin', 'GC', 'Length', 'Depth']],
                                        on='Original_Contig_Name')

#df2=pd.merge(df,df_scaff, right_on='Original_Contig_Name',
#	left_index=True, how='left') 




# Create output file  with the mapping file 


outfilename = args.filename + "_mapping_scaffolds.tsv"


df2.to_csv(outfilename, sep ="\t",header=True)


# Create output file with total numbers 

outfile_stats  = args.filename + "_stats.tsv"

df3=df2.nunique()

df3.to_csv(outfile_stats, sep="\t",header=True) 

#Create a file only with matching databases 

outfile_taxonomy  = args.filename + "_dropna.tsv"

taxonomy=df2[['Original_Contig_Name', 'Gene_Type','Gene_Length','Lineage','COG_ID','Cog_%ID','PFAM_ID','KO_Term','Bin']]
taxonomy=taxonomy.dropna()
taxonomy.to_csv(outfile_taxonomy, sep ="\t",header=True)


outfile_kos  = args.filename + "_kosperbin.tsv"


#create a file with bins and kos per bin 
group = taxonomy.groupby('Bin')
kos = pd.DataFrame(group.apply(lambda x: x['KO_Term'].unique()))
kos.to_csv(outfile_kos,sep="\t")



outfile_cogs  = args.filename + "_cogsperbin.tsv"


#create a file with bins and cogs per bin
group = taxonomy.groupby('Bin')
cogs = pd.DataFrame(group.apply(lambda x: x['COG_ID'].unique()))
cogs.to_csv(outfile_cogs,sep="\t")


outfile_pfams  = args.filename + "_pfamsperbin.tsv"
#create a file with bins and pfams per bin

group = taxonomy.groupby('Bin')
pfams = pd.DataFrame(group.apply(lambda x: x['PFAM_ID'].unique()))
pfams.to_csv(outfile_pfams,sep="\t")


