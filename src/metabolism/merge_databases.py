#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     merge_databases.py
# Purpose:  merge the metabolic databases together
# @uthors:      sbs - sbs2756@utexas.edu
# Created:     2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------

import argparse
import glob
import os
import pandas as pd


def main():

	#Command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--databases_path', required=True, help=("Path for directory containing all the database files."))
	parser.add_argument('-d', '--database', required=True, help=("Database(s) of interest to merge together. Input as a list."))
	args = parser.parse_args()

	#Create output directory if not already present
	if not os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/../../output"):
		os.makedirs(os.path.dirname(os.path.abspath(__file__)) + "/../../output")

	dfs_list = []
	files_list = []

	#Make sure directory path is valid
	if args.databases_path[-1] != '/':
		databases_path = args.databases_path + '/'
	else:
		databases_path = args.databases_path

	#Verify database name is valid
	if ' ' in args.database:
		databases_list = args.database.upper().strip('[]').split(', ')
	else:
		databases_list = args.database.upper().strip('[]').split(',')
	for check in databases_list:
		if check not in ['KEGG', 'COG', 'PFAM', 'EC_NUMBER']:
			print('Invalid database name was entered!')
			quit()

	#Get all dataframes and merge them
	for d in databases_list:
		merged_list = []
		for file in glob.glob(databases_path + d + '_metabolic_profile*'):
			database_df = pd.read_csv(file, sep='\t', index_col=False)
			merged_list.append(database_df)
		
		#Save each merged dataframe
		merged_dfs = pd.concat(merged_list)
		merged_dfs = merged_dfs.groupby(by=d).sum()
		merged_dfs = merged_dfs.reset_index()
		files_list.append('merged_' + d + '_metabolic_profile.tsv')
		merged_dfs.to_csv(os.path.dirname(os.path.abspath(__file__)) + '/../../output/merged_' + d + '_metabolic_profile.tsv', sep='\t', index=False)
		if len(databases_list) > 1:
			dfs_list.append(merged_dfs)

	if len(databases_list) > 1:
		#Combine multiple database dataframes into one
		final_df = pd.concat(dfs_list)
		final_df = final_df.set_index(databases_list)
		final_df = final_df.fillna(0.0)
		#Save dataframe containing multiple database columns
		files_list.append('final_merged_metabolic_profile.tsv')
		final_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + '/../../output/final_merged_metabolic_profile.tsv', sep='\t', index=True)

	print("Success!\nThe following files have been saved in the \"output\" directory:\n")
	for f in files_list:
		print(f)
	print()


if __name__ == '__main__':
	main()