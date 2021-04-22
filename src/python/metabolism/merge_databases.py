# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     merge_databases.py
# Purpose:  merge the metabolic databases together
#
# @uthors:      sbs - sbs2756@utexas.edu
#
# Created:     2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------

import argparse
import glob
import pandas as pd

def main():

	#Command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--databases_path', required=True, help=("Path for directory containing all the database files."))
	parser.add_argument('-d', '--database', required=True, help=("Database(s) of interest to merge together. Input as a list."))
	args = parser.parse_args()

	dfs_list = []

	#Make sure directory path is valid
	if args.databases_path[-1] != '/':
		databases_path = args.databases_path + '/'
	else:
		databases_path = args.databases_path

	#Verify database name is valid
	databases_list = args.database.upper().strip('[]').split(', ')
	for check in databases_list:
		if check not in ['KEGG', 'COG', 'PFAM', 'EC_NUMBER']:
			print('Invalid database name was entered!')
			quit()

	#Get all dataframes and merge them
	for d in databases_list:
		counter = 1
		for file in glob.glob(databases_path + d + '_metabolic_profile*'):
			if counter == 1:
				merged_dfs = pd.read_csv(file, sep='\t', index_col=False)
			else:
				f = pd.read_csv(file, sep='\t', index_col=False)
				merged_dfs = merged_dfs.merge(f, on=d, how='outer')
			counter+=1

		#Save each merged dataframe
		merged_dfs.to_csv('../../../output/merged_' + d + '_metabolic_profile.tsv', sep='\t', index=False)
		if len(databases_list) > 1:
			dfs_list.append(merged_dfs)

	if len(databases_list) > 1:
		#Combine multiple database dataframes into one
		final_df = pd.concat(dfs_list)
		final_df = final_df.set_index(databases_list)
		#Save dataframe containing multiple database columns
		final_df.to_csv('../../../output/final_merged_metabolic_profile.tsv', sep='\t', index=True)


if __name__ == '__main__':
	main()