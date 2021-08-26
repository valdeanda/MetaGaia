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
import glob
import numpy as np
import os
import pandas as pd
import re

def format_vibrant(vibrant_df, vibrant_name):
	"""
	This function formats the Vibrant scaffold names by adding the sampling_site name as a prefix.
	Input(s):
	vibrant_df is a pandas dataframe of the Vibrant output.
	vibrant_name is a string of the name of the vibrant file.
	Output(s):
	vibrant_df is a pandas datafrmae that has been modified to include a Sample column and an altered scaffold column.
	"""

	sample_name = re.search('individuals_(.+?)_scaffold', vibrant_name).group(1)
	vibrant_df['Sample'] = sample_name
	vibrant_df['scaffold'] = sample_name + '_' + vibrant_df['scaffold']

	return vibrant_df


def main():

	#Command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-p','--pathway_database', required=True, help="Input file in tsv format that contains the counts of each pathway per bin.")
	parser.add_argument('-vf', '--vibrant_file', required=False, help="Input file from Vibrant in tsv format that contains phage scaffolds mapped to metabolic pathways.")
	parser.add_argument('-vp', '--vibrant_path', required=False, help="Input path from Vibrant files in tsv format that contains phage scaffolds mapped to metabolic pathways.")
	parser.add_argument('-d', '--database', required=True, help="Input only one of the following database names the user is interested in analyzing: KEGG or PFAM.")
	args = parser.parse_args()

	#Create output directory if not already present
	if not os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/../../output"):
		os.makedirs(os.path.dirname(os.path.abspath(__file__)) + "/../../output")

	drop_list = []

	print("Beginning to analyze phage and host metabolisms.")

	#Read in input files
	database = args.database.upper()
	host_df = pd.read_csv(args.pathway_database, sep='\t', index_col=False)
	#Drop unneeded database columns
	if 'EC_NUMBER' in host_df.columns.tolist():
		host_df = host_df.drop(columns=['EC_NUMBER'])
	if 'COG' in host_df.columns.tolist():
		host_df = host_df.drop(columns=['COG'])
	if 'PFAM' in host_df.columns.tolist() and database != 'PFAM':
		host_df = host_df.drop(columns=['PFAM'])
	elif 'KEGG' in host_df.columns.tolist() and database != 'KEGG':
		host_df = host_df.drop(columns=['KEGG'])
	host_df = host_df.drop(columns=['NoBin', 'Total'])
	host_df = host_df.dropna()

	#Read in vibrant output
	if args.vibrant_file:
		phage_df = pd.read_csv(args.vibrant_file, sep='\t')
		phage_df = format_vibrant(phage_df, args.vibrant_file.split('/')[-1])
	else:
		vibrant_list = []

		if args.vibrant_path[-1] != '/':
			vibrant_path = args.vibrant_path + '/'
		else:
			vibrant_path = args.vibrant_path
		#Format each vibrant file in directory and concatenate
		for vibrant in glob.glob(vibrant_path + "VIBRANT_AMG*"):
			vibrant_df = pd.read_csv(vibrant, sep='\t', index_col=False)
			vibrant_df = format_vibrant(vibrant_df, vibrant.split('/')[-1])
			vibrant_list.append(vibrant_df)
		phage_df = pd.concat(vibrant_list)

	#Format metabolic profile into two columns: KEGG pathways and Bins
	host_df = pd.melt(host_df, id_vars=database, value_vars=host_df.columns.tolist()[1:])
	host_df = host_df.rename(columns={'variable': 'Bin'})
	#Drop rows that do not have KEGG pathways in bins
	host_df = host_df[host_df.value != 0]
	host_df = host_df.drop(columns=['value'])

	#Add Sample column to metabolic profile (host dataframe)
	host_df['Sample'] = host_df['Bin'].str.extract('(.+?)_Bin', expand=False)

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
	phage_host_df = host_df.merge(phage_df[['scaffold', 'Sample', database]], on=['Sample', database], how='outer')
	phage_host_df = phage_host_df.fillna(np.nan)

	print('Determining if phages and hosts have common metabolic pathways. This may take a while. Ignore any warnings.')
	phage_host_df = phage_host_df.drop(columns=['Sample'])
	#Create new column determining presence of metabolic pathway
	phage_host_df['Presence'] = np.nan
	#Subset dataframe to alter "Presence" column accordingly
	#Scaffold dataframe
	dfp = phage_host_df[(phage_host_df['Bin'].isnull()) & (phage_host_df['scaffold'].notnull())]
	dfp['Presence'] = 'phage'
	#Bins dataframe
	dfh = phage_host_df[(phage_host_df['scaffold'].isnull()) & (phage_host_df['Bin'].notnull())]
	dfh['Presence'] = 'host'
	#Dataframe containing both bins and scaffolds
	dfb = phage_host_df[(phage_host_df['scaffold'].notnull()) & (phage_host_df['Bin'].notnull())]
	dfb['Presence'] = 'both'
	phage_host_df = pd.concat([dfp, dfh, dfb])

	#Drop all rows with NaN in the Presence column
	phage_host_df = phage_host_df.dropna(subset=['Presence'])
	#Save file
	phage_host_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + '/../../output/phage_host_metabolism.tsv', sep='\t', index=False)

	#Only keep the bins mapped to scaffolds
	if database == 'KEGG':
		phage_host_df = phage_host_df[phage_host_df['Presence'] == 'both'].drop(columns=['Presence', 'KEGG'])
	elif database == 'PFAM':
		phage_host_df = phage_host_df[phage_host_df['Presence'] == 'both'].drop(columns=['Presence', 'PFAM'])
	#Save file
	phage_host_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + '/../../output/phage_host_mapping.tsv', sep='\t', index=False)
	print("Success!\nThe following files have been saved in the \"output\" directory:\n\nphage_host_metabolism.tsv\nphage_host_mapping.tsv\n")

if __name__ == "__main__":
	main()


