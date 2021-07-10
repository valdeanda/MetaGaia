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
import copy
import glob
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
		self.parser.add_argument('-i', '--imganno_file', required=False, help=("Input file or files in tsv format. Rows are genes columns are IMG annotations."))
		self.parser.add_argument('-p', '--imganno_path', required=False, help=("Input path containing only the IMG files. Rows are genes columns are IMG annotations in each file."))
		self.parser.add_argument('-m', '--mapping', required=True, help=("A tsv file containing original contig name, sample, and bin columns. Created from the files_prep.py script."))
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
	files_list is a list containing all the files created in this function.
	"""

	files_list = []

	#Map scaffolds 
	mapping_df = mapping_df[['Original_Contig_Name', 'Bin']]
	if 'Bin' not in img_df.columns:
		scaff_df = pd.merge(img_df, mapping_df, on='Original_Contig_Name', how='left')
	else:
		scaff_df = copy.deepcopy(img_df)
	#Rename columns
	scaff_df = scaff_df.rename(columns={'Original_Contig_Name': 'Scaffold'})
	mapping_df = mapping_df.rename(columns={'Original_Contig_Name': 'Scaffold'})

	#Extract columns of interest
	taxonomy_df = scaff_df[['Scaffold','COG_ID','PFAM_ID','KO_Term', 'EC_Number', 'Bin']]

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
	ec_df = pd.DataFrame(taxonomy_df.apply(lambda x: x['EC_Number'].unique().tolist()), columns=['EC_NUMBER'])
	#Merge KEGG/COG/PFAM with EC dataframe
	databases_df = databases_df.merge(ec_df, on='Scaffold', how='outer')
	#Map bins to each scaffold
	databases_df = databases_df.merge(mapping_df, on='Scaffold', how='left')
	#Fill bin NaNs as NoBin
	databases_df = databases_df.fillna({'Bin': 'NoBin'})

	files_list.append(uniquify(os.path.dirname(os.path.abspath(__file__)) + "/../../output/mapped_scaffolds.tsv").split('/')[-1])
	#Save file in output folder
	databases_df.to_csv(uniquify(os.path.dirname(os.path.abspath(__file__)) + "/../../output/mapped_scaffolds.tsv"), index=False, sep="\t")

	return [databases_df, files_list]


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
		path = filename + str(i) + extension
		i+=1

	return path


def get_database_counts(extract_list, databases_df):
	"""
	This function saves a file containing the number of times a database value appears in a bin.
	Input(s):
	extract_list is a list of the desired databases to analyze.
	databases_df is a pandas dataframe that maps each scaffold to metabolic pathway that are found for it.
	Output(s):
	files_list is a list containing all the files created in this function.
	"""

	files_list = []
	dfs_list = []

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
		#Rename index
		count_df.index.names = [d]
		count_df = count_df.reset_index()

		dfs_list.append(count_df)

	#Combine multiple database dataframes into one
	final_df = pd.concat(dfs_list)
	final_df = final_df.set_index(extract_list)
	final_df = final_df.fillna(0.0)
	#Save dataframe containing multiple database columns
	files_list.append(uniquify(os.path.dirname(os.path.abspath(__file__)) + '/../../output/complete_metabolic_profile.tsv').split('/')[-1])
	final_df.to_csv(uniquify(os.path.dirname(os.path.abspath(__file__)) + '/../../output/complete_metabolic_profile.tsv'), sep='\t', index=True)

	return files_list

def concat_img(img_path):
	"""
	This function will concatenate the IMG files if there are more than one present within a directory.
	Input(s):
	img_path is a string containing the path to the img files needed for concatenation.
	Output(s):
	img_file is a pandas dataframe that contains all the IMG information in one file.
	"""

	files_list = []
	img_list = []

	if img_path[-1] != '/':
		img_path = img_path + '/'

	for img in glob.glob(img_path + '*'):
		if img.endswith('.tsv'):
			img_df = pd.read_csv(img, sep='\t')
		elif img.endswith('.csv'):
			img_df = pd.read_csv(img)
		elif img.endswith('.txt'):
			img_df = pd.read_csv(img, sep='\s+')
		img_list.append(img_df)

	img_file = pd.concat(img_list)
	files_list.append(uniquify(os.path.dirname(os.path.abspath(__file__)) + '/../../output/IMG_consolidated_master.tsv').split('/')[-1])
	img_file.to_csv(uniquify(os.path.dirname(os.path.abspath(__file__)) + '/../../output/IMG_consolidated_master.tsv'), sep='\t')

	return [img_file, files_list]


def main():

	#Command line arguments
	arguments = Command_line_args()

	saved_files = []

	#Create output directory if not already present
	if not os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/../../output"):
		os.makedirs(os.path.dirname(os.path.abspath(__file__)) + "/../../output")

	#If the database the user wants is not present, quit the program
	if ' ' in arguments.args.database:
		databases_list = arguments.args.database.upper().strip('[]').split(', ')
	else:
		databases_list = arguments.args.database.upper().strip('[]').split(',')
	for check in databases_list:
		if check not in ['KEGG', 'COG', 'PFAM', 'EC_NUMBER']:
			print('The ' + check + ' column is not mapped in the mapped_scaffolds file!')
			quit()

	#Read mapping file created previously (bin_abundance step)
	mapping_df = pd.read_csv(arguments.args.mapping, sep="\t")

	print('Creating a consolidated IMG file containing all the IMG information for ALL sample(s).')
	#Read in IMG annotated file
	if arguments.args.imganno_file:
		if arguments.args.imganno_file.endswith('.tsv'):
			img_df = pd.read_csv(arguments.args.imganno_file, sep='\t')
		elif arguments.args.imganno_file.endswith('.csv'):
			img_df = pd.read_csv(arguments.args.imganno_file)
		elif arguments.args.imganno_file.endswith('.txt'):
			img_df = pd.read_csv(arguments.args.imganno_file, sep='\s+')
	else:
		img_concat = concat_img(arguments.args.imganno_path)
		img_df = img_concat[0]
		saved_files = img_concat[1]
	print('Beginning to map scaffolds to bins and database values. Getting count of each database value in each bin. This may take a while.')
	#Map scaffolds to each of the databases
	scaffold_mapping = map_scaffolds(arguments, img_df, mapping_df)
	databases_df = scaffold_mapping[0]
	if len(saved_files) == 0:
		saved_files = scaffold_mapping[1]
	else:
		saved_files = saved_files + scaffold_mapping[1]
	#Get the counts of each metabolic pathway in each bin
	saved_files = saved_files + get_database_counts(databases_list, databases_df)
	
	print("Success!\nThe following files have been saved in the \"output\" directory:\n")
	for f in saved_files:
		print(f)
	print()
	

if __name__ == '__main__':
	main()




	






