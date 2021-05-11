#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     pathway_extraction.py
# Purpose:  use the mapped scaffolds to extract the pathways of interest
#
# @uthors:     sbs - sbs2756@utexas.edu
#
# Created:     2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------


import argparse
import os
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))


class Command_line_args():
	"""
    This class contains all the arguments the user inputs for the class to run.
    Input(s):
    No other inputs needed.
    Output(s):
    None.
    """

	def __init__(self):

		#Command line arguments
		self.parser = argparse.ArgumentParser()
		self.parser.add_argument('-c', '--customdata', required=True, help=("Input file in tsv format with at least one column with the header containing the pathways of interest."))
		self.parser.add_argument('-p', '--pathway_database', required=True, help=("Input file in tsv format that contains the counts of each pathways per bin."))
		self.parser.add_argument('-d', '--database', required=True, help=("Database(s) of interest to merge together. Input as a list."))
		self.args = self.parser.parse_args()


def main():

	print('Make sure the pathway codes in the user provided file are formatted the same as it is throughout MetaGaia!!')
	#Command line arguments
	arguments = Command_line_args()
	pathways_lst = []

	#Read in files
	user_df = pd.read_csv(arguments.args.customdata, sep='\t', index_col=False)
	pathway_df = pd.read_csv(arguments.args.pathway_database, sep='\t', index_col=False)

	#Verify database name is valid
	databases_list = args.database.upper().strip('[]').split(', ')
	for check in databases_list:
		if check not in ['KEGG', 'COG', 'PFAM', 'EC_NUMBER']:
			print('Invalid database name was entered!')
			quit()
	else:
		#Convert to list
		user_df[check] = user_df[check].str.split(', ')
		index_cols = [col for col in user_df.columns if col != check]
		user_df = (user_df.set_index(index_cols)[check].apply(pd.Series).stack().reset_index().drop('level_'+str(len(index_cols)), axis=1).rename(columns={0:check}))
		for val in user_df[check]:
			if check == 'KEGG':
				if val[:3] != 'KO:':
					user_df.loc[user_df[check] == val, check] = 'KO:' + val
			elif check == 'PFAM':
				if val[:4] != 'pfam':
					user_df.loc[user_df[check] == val, check] = 'pfam' + val[2:]
			elif check == 'EC_NUMBER':
				if val[:3] != 'EC:':
					user_df.loc[user_df[check] == val, check] = 'EC:' + val

	print('Beginning to extract pathway information.')
	for d in databases_list:
		#Subset dataframe
		extracted_df = pathway_df[pathway_df[d].isin(user_df[d])]
		#If user provided file has multiple columns, map them to the extracted pathway dataframe
		remove_cols = [col for col in databases_list if col != d]
		extracted_df = extracted_df.drop(columns=remove_cols)
		if len(user_df.drop(columns=databases_list).columns.tolist()) > 1:
			extracted_df = extracted_df.merge(user_df.drop(columns=remove_cols), on=d, how='left')
		#Rename datbase column so that concatenation can occur
		cols_list = extracted_df.columns.tolist()
		cols_list[cols_list.index(d)] = 'Database'
		extracted_df.columns = cols_list
		#Add dataframes to list
		pathways_lst.append(extracted_df)

	#Combine all dataframes
	final_df = pd.concat(pathways_lst)
	final_df = final_df.set_index('Database')
	#Save file
	final_df.to_csv('../../output/extracted_pathways.tsv', sep='\t', index=True)
	print('Finished extracting information! File is saved in the \"output\" folder as \"extracted_pathways.tsv\".')


if __name__ == '__main__':
	main()





