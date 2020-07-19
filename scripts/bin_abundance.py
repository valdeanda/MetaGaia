
#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     Bin_abundance.py
# Purpose:  calculate the abundance of bins based of scaffold depth
#
# @uthors:      vda - valdeanda@utexas.edu
#
# Created:     2020
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------


import argparse
import pandas as pd


# options
epilog = """Example:

$ python3 bin_abundance.py  -l lenghtpersample.tab  -m GB_CP_MappingFile4.tab  -t Taxonomy.tsv"""

parser = argparse.ArgumentParser(description=__doc__, epilog=epilog)
parser.add_argument('-l', '--length', required=True,
                    help=("Total lenght of the samples. Rows are samples "
                    "column is length ."))

parser.add_argument('-m', '--mappingfile', required=True,
                    help=("Tabular file containing the following columns "
                          "Original_Contig_Name    Bin GC  Length  Depth   Sample "
                          "Make sure that you have those exact columns"))


parser.add_argument('-t', '--taxonomy', required=True,
                    help=("Tabular file with Bins taxonomy "
                   "Bin, Proposed_Name, Taxonomy, Superphylum and Confirmed Taxonomy "
		   "are columns"))



args = parser.parse_args()


length       = args.length
mapping_file = args.mappingfile
taxonomy     = args.taxonomy

df_length   = pd.read_csv(length,sep="\t")
df_mapping  = pd.read_csv(mapping_file,sep="\t",index_col=False)
df_taxonomy = pd.read_csv(taxonomy, sep ="\t")


# Map the total lengt in the mapping file

df2=pd.merge(df_mapping,df_length[['Sample', 'total_length(bp)']], on='Sample')

# Calculate coverage

df2['Cov']= df2['Length'] * df2['Depth']/df2['Length']

# Calculate abundance
df2['Abund']=df2['Cov']/df2['total_length(bp)']
#add taxonomy to the mapping file
#df3=pd.merge(df2,df_taxonomy[['Bin', 'Proposed_Name','Taxonomy','Superphylum','Confirmed Taxonomy']], on='Bin')
df3=pd.merge(df2,df_taxonomy[['Bin', 'Domain','Class','Proposed taxonomy','Superphylum/Phylum','Notes']], on='Bin')
#abundance of each bin

df3['Abun_bin'] = df3['Abund'].groupby(df3['Bin']).transform('sum')

outtaxonomy = args.taxonomy + "_abundance.tsv"
df_final=df3.set_index('Original_Contig_Name')
df3.to_csv(outtaxonomy,sep='\t')


print("[END] Done computing abundance....................\n"
      "Please check the output file:\n",
      "Mapping File with Bin Abundance:", args.taxonomy + "_abundance.tsv\n")
