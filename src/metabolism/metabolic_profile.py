#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     metabolic_profile.py
# Purpose:  map scaffolds to their respective bins and determine the number of times a metabolic pathway appears in within each bin
# @uthors:      sbs - sbs2756@utexas.edu
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
		self.parser.add_argument('-i', '--imganno', required=True, help=("Input file in tsv format. Rows are genes columns are IMG annotations."))
		self.parser.add_argument('-m', '--mapping', required=True, help=('A tsv file containing original contig name, sample, and bin columns'))
		self.parser.add_argument('-d', '--database', required=True, help=("Database(s) of interest to merge together. Input as a list."))
		self.parser.add_argument('-c', '--consistency', required=False, default="", help=('Boolean value that determines if scaffolds not containing a value for each database should be kept. Leave blank if consistency check is not needed.'))
		self.args = self.parser.parse_args()


def map_scaffolds(arguments, img_df, mapping_df):
	"""
	This function saves and returns a dataframe that contains the KEGG, COG, PFAM, and EC_Number values for each scaffold.
	Input(s):
	arguments is a class containing all the command line arguments.
	img_df is a pandas dataframe containing the IMG annotations.
	mapping_df is a pandas dataframe mapping the original contig name to its corresponding bin and sample.
	Output(s):
	databases_df is a pandas dataframe that maps each scaffold to metabolic pathways that are found for it.
	"""

	#Map scaffolds 
	mapping_df = mapping_df[['Original_Contig_Name', 'Bin']]
	scaff_df = pd.merge(img_df, mapping_df, on='Original_Contig_Name')
	#Rename columns
	mapping_df.columns = ['Scaffold', 'Bin']

	#Extract columns of interest
	taxonomy_df = scaff_df[['Original_Contig_Name','COG_ID','PFAM_ID','KO_Term', 'EC_Number', 'Bin']]
	taxonomy_df.columns = ['Scaffold', 'COG_ID', 'PFAM_ID', 'KO_Term', 'EC_NUMBER', 'Bin']

	#Drop row if NaN is present in any of the four databases
	if arguments.args.consistency:
		taxonomy_df = taxonomy_df.dropna()


	#Create a dataframe with bins and kos per bin 
	taxonomy_df = taxonomy_df.groupby('Scaffold')
	kegg_df = pd.DataFrame(taxonomy_df.apply(lambda x: x['KO_Term'].unique().tolist()), columns=['KEGG'])


	#Create a dataframe with bins and cogs per bin
	cog_df = pd.DataFrame(taxonomy_df.apply(lambda x: x['COG_ID'].unique().tolist()), columns=['COG'])
	#Merge KEGG and COG dataframes
	databases_df = kegg_df.merge(cog_df, on='Scaffold', how='outer')

	#Create a dataframe with bins and pfams per bin
	pfam_df = pd.DataFrame(taxonomy_df.apply(lambda x: x['PFAM_ID'].unique().tolist()), columns=['PFAM'])
	#Merge KEGG/COG with PFAM dataframe
	databases_df = databases_df.merge(pfam_df, on='Scaffold', how='outer')

	#Create a dataframe with bins and ECs per bin
	ec_df = pd.DataFrame(taxonomy_df.apply(lambda x: x['EC_NUMBER'].unique().tolist()), columns=['EC_NUMBER'])
	#Merge KEGG/COG/PFAM with EC dataframe
	databases_df = databases_df.merge(ec_df, on='Scaffold', how='outer')
	#Map bins to each scaffold
	databases_df = databases_df.merge(mapping_df, on='Scaffold', how='left')

	#Save file in output folder
	databases_df.to_csv(uniquify("../../output/mapped_scaffolds.tsv"), index=False, sep="\t")

	return databases_df


def uniquify(path):
	"""
	This function validates that the file is a unique filename in the directory.
	Input(s):
	path is a string containing the path to the file.
	Output(s):
	path is a string containing the path to a unique filename.
	"""

	filename, extension = os.path.splitext(path)
	i = 1

	while os.path.exists(path):
		path = filename + str(counter) + extension
		i+=1

	return path


def get_database_counts(extract_list, databases_df):
	"""
	This function saves a file containing the number of times a database value appears in a bin.
	Input(s):
	extract_list is a list of the desired databases to analyze.
	databases_df is a pandas dataframe that maps each scaffold to metabolic pathway that are found for it.
	Output(s):
	count_df is a pandas dataframe that tracks the count of each metabolic pathway per bin.
	"""

	#For each database the user wants to analyze
	for d in extract_list:
		db_list = []
		bin_list = []
		#For the column containing the desired database
		for row in range(len(databases_df[d])):
			#For the row in the database column
			for k in range(len(databases_df[d][row])):
				#Append metabolic pathway and bin in a list of values per scaffold
				db_list.append(databases_df[d][row][k])
				bin_list.append(databases_df['Bin'][row])

		#Create dataframe with each database value mapped to the bin that contains it
		db_bin_df = pd.DataFrame(data={d: db_list, 'Bin': bin_list})
		#Group by and count the number of times the metabolic pathway appears in a bin
		count_df = db_bin_df.groupby([d, 'Bin']).size().reset_index()
		#Rename columns
		count_df.columns = [d, 'Bin', 'Count']
		#Pivot long to wide so that row index are metabolic pathway and columns are the bin names
		count_df = count_df.pivot(index=d, columns='Bin')
		#Reorganize column names
		count_df.columns = count_df.columns.droplevel()
		#Replace NaN with 0.0
		count_df = count_df.fillna(0)
		#Create a total columns
		count_df['Total'] = count_df.sum(axis=1)

		#Save file in the output folder
		count_df.to_csv(uniquify('../../output/' + d + '_metabolic_profile.tsv'), index=True, sep='\t')


def main():

	#Command line arguments
	arguments = Command_line_args()
	#Read in IMG annotated file
	img_df = pd.read_csv(arguments.args.imganno, sep="\t")
	#Read mapping file created previously (bin_abundance step)
	mapping_df = pd.read_csv(arguments.args.mapping, sep="\t")

	print('Beginning to map scaffolds to bins and database values.')
	#Map scaffolds to each of the databases
	databases_df = map_scaffolds(arguments, img_df, mapping_df)
	print('Finished mapping scaffolds!')

	#If the database the user wants is not present, quit the program
	databases_list = args.database.upper().strip('[]').split(', ')
	for check in databases_list:
		if check not in ['KEGG', 'COG', 'PFAM', 'EC_NUMBER']:
			print('The requested database is not mapped in the mapped_scaffolds file!')
			quit()

	print('Getting count of each database value in each bin.')
	#Get the counts of each metabolic pathway in each bin
	get_database_counts(databases_list, databases_df)
	print('Finished getting counts! All outputs should be saved in the \"output\" folder.')
	

if __name__ == '__main__':
	main()




	






