#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     id_common_metabolism.py
# Purpose:  identify common metabolic pathways present within phages (scaffolds) and their hosts (bins)
# @uthors:     sbs - sbs2756@utexas.edu
# Created:     April 2020
# Last updated April  2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------

import argparse
import numpy as np
import os
import pandas as pd


def main():

	#Command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--pathway_database', required=True, help="Input file in tsv format that contains the counts of each pathways per bin.")
	parser.add_argument('-v', '--vibrant_outout', required=True, help="Input file from Vibrant in tsv format that contains phage scaffolds mapped to metabolic pathways.")
	parser.add_argument('-d', '--database', required=True, help="Input only one of the following database names the user is interested in analyzing: KEGG, COG, PFAM, or EC_NUMBER.")
	args = parser.parse_args()

	drop_list = []

	print("Beginning to analyze phage and host metabolisms.")

	#Read in input files
	host_df = pd.read_csv(args.pathway_database, sep='\t')
	host_df = host_df.drop(columns=['NoBin'])
	phage_df = pd.read_csv(args.vibrant_outout, sep='\t')
	database = args.database.upper()

	host_df = pd.melt(host_df, id_vars=database, value_vars=host.columns.tolist()[1:])
	host_df = host_df.rename(columns={'variable': 'Bin'}).drop(columns=['value'])

	#Rename database columns
	if database == 'KEGG':
		phage_df = phage_df.rename(columns={'AMG KO': 'KEGG'})
	elif database == 'Pfam':
		phage_df = phage_df.rename(columns={'Pfam': 'PFAM'})
	#Change pathway name format
	for val in phage_df[database]:
		if database == 'KEGG':
			if val[:3] != 'KO:':
				phage_df.loc[phage_df[database] == val, database] = 'KO:' + val
		elif database == 'PFAM':
			if val[:4] != 'pfam':
				phage_df.loc[phage_df[database] == val, database] = 'pfam' + val[2:]

	#Join the two dataframes
	phage_host_df = host_df.merge(phage_df[['scaffold', database]], on=database, how='outer')
	phage_host_df = phage_host_df.fillna('None')

	print('Determining if phages and hosts have common metabolic pathways.')
	#Create new column determining presence of metabolic pathway
	phage_host_df['Presence'] = np.nan
	for i in range(len(phage_host_df)):
		if phage_host_df['Bin'].iloc[i] == 'None' and hage_host_df['scaffold'].iloc[i] == 'None':
			continue
		elif phage_host_df['Bin'].iloc[i] == 'None':
			phage_host_df['Presence'].iloc[i] = 'phage'
		elif phage_host_df['scaffold'].iloc[i] == 'None':
			phage_host_df['Presence'].iloc[i] = 'host'
		else:
			phage_host_df['Presence'].iloc[i] = 'both'
	#Drop all rows with NaN in the Presence column
	phage_host_df = phage_host_df.dropna()
	#Save file
	phage_host_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + '/../../output/phage_host_metabolism.tsv', sep='\t', index=False)

	#Only keep the bins mapped to scaffolds
	phage_host_df = phage_host_df[phage_host_df['Presence'] == 'both'].drop(columns=['Presence', 'KEGG'])
	#Save file
	phage_host_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + '/../../output/phage_host_mapping.tsv', sep='\t', index=False)
	print('Success! All outputs are saved in the \"output\" directory.')

if __name__ == "__main__":
	main()


