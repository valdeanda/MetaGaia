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
		self.parser.add_argument('-c', '--customdata', required=True, help=("Input file in tsv format with at least one column with the header containing the pathways of interest. The headers of the metabolic pathway must correspond with its respective database (eg. KEGG, COG, PFAM, or EC_NUMBER)."))
		self.parser.add_argument('-p', '--pathway_database', required=True, help=("Input file in tsv format that contains the counts of each pathways per bin."))
		self.parser.add_argument('-d', '--database', required=True, help=("Database(s) of interest to merge together. Input as a list."))
		self.args = self.parser.parse_args()


def main():

	print('Make sure the pathway codes in the user provided file are formatted the same as it is throughout MetaGaia!')
	#Command line arguments
	arguments = Command_line_args()

	#Create output directory if not already present
	if not os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/../../output"):
		os.makedirs(os.path.dirname(os.path.abspath(__file__)) + "/../../output")

	pathways_lst = []

	#Read in files
	user_df = pd.read_csv(arguments.args.customdata, sep='\t', index_col=False)
	pathway_df = pd.read_csv(arguments.args.pathway_database, sep='\t', index_col=False)

	#Verify database name is valid
	if ' ' in arguments.args.database:
		databases_list = arguments.args.database.upper().strip('[]').split(', ')
	else:
		databases_list = arguments.args.database.upper().strip('[]').split(',')
	for check in databases_list:
		if check not in ['KEGG', 'COG', 'PFAM', 'EC_NUMBER']:
			print('Invalid database name was entered!')
			quit()
		else:
			#Unstack columns with muliple values if necessary
			if user_df[check].str.contains(',').sum() > 1:
				if user_df[check].str.contains(', ').sum() > 1:
					add_df = user_df[check].str.split(', ').apply(pd.Series,1).stack()
				else:
					add_df = user_df[check].str.split(',').apply(pd.Series,1).stack()
				add_df.index = add_df.index.droplevel(-1)
				add_df.name = check
				del user_df[check]
				user_df = user_df.join(add_df)

			for val in user_df[check]:
				if check == 'KEGG':
					if str(val)[:3] != 'KO:':
						user_df.loc[user_df[check] == str(val), check] = 'KO:' + str(val)
				elif check == 'PFAM':
					if str(val)[:4] != 'pfam':
						user_df.loc[user_df[check] == str(val), check] = 'pfam' + str(val)[2:]
				elif check == 'EC_NUMBER':
					if str(val)[:3] != 'EC:':
						user_df.loc[user_df[check] == str(val), check] = 'EC:' + str(val)
	user_df = user_df.reset_index().drop(columns=['index'])

	print('Beginning to extract pathway information.')
	for d in databases_list:
		#Subset dataframe
		extracted_df = pathway_df[pathway_df[d].isin(user_df[d].unique())]
		#If user provided file has multiple columns, map them to the extracted pathway dataframe
		remove_cols = [col for col in databases_list if col != d]
		extracted_df = extracted_df.drop(columns=remove_cols)
		if len(user_df.drop(columns=databases_list).columns.tolist()) > 1:
			extracted_df = extracted_df.merge(user_df.drop(columns=remove_cols), on=d, how='left')
		#Rename datbase column so that concatenation can occur
		cols_list = extracted_df.columns.tolist()
		cols_list[cols_list.index(d)] = 'Database'
		extracted_df.columns = cols_list
		extracted_df = extracted_df.dropna().drop_duplicates()
		#Add dataframes to list
		pathways_lst.append(extracted_df)

	#Combine all dataframes
	final_df = pd.concat(pathways_lst)
	final_df = final_df.set_index('Database')
	#Save file
	final_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + '/../../output/extracted_pathways.tsv', sep='\t', index=True)
	print("Success!\nThe following files have been saved in the \"output\" directory:\n\nextracted_pathways.tsv\n")


if __name__ == '__main__':
	main()





