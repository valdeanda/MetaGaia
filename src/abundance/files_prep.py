#!/usr/bin/env python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:        files_prep.py
# Purpose:     create the input files needed to calculate bin abundacnces.
# @uthors:     sbs - sahil2699@gmail.com
# Created:     2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
import argparse
import copy
import glob
import numpy as np
import os
import pandas as pd
import re


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
		self.parser.add_argument('-b', '--bin_samples', required=True, type=str, help='Input file containing 2 columns: 1st column contains bin names and 2nd column contains sample names. Requires header for each column, \"Bin\" and \"Sampling_Site\" respectively. File must have tsv, csv, or txt file extension.')
		self.parser.add_argument('-d', '--depth', required=False, default="", type=str, help='Input file containing 5 columns: contigName, contigLen, totalAvgDepth, Sample, and Depth. This file can be created with jgi_summarize_bam_contig_depths. File must have tsv, csv, or txt file extension.')
		self.parser.add_argument('-f', '--depth_dir', required=False, default="", type=str, help='Input directory path containing all the depth txt files that must be concatenated. No other txt files should be present within the directory!')
		self.parser.add_argument('-q', '--fastq_dir', required=True, type=str, help='Input directory path containing all the fastq files.')
		self.parser.add_argument('-a', '--fna_dir', required=True, type=str, help='Input directory path containing all the fna files.')
		self.args = self.parser.parse_args()


def create_reads_file(arguments):
	"""
	This function creates the reads file needed.
	Input(s):
	arguments are the command line arguments.
	Output(s):
	A file containing the reads information is saved in the "output" folder.
	"""

	if os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/../../output/fastq_read_counts.txt"):
		os.remove(os.path.dirname(os.path.abspath(__file__)) + "/../../output/fastq_read_counts.txt")
	#Run bash script to create a text file containing the reads information
	bash_command = "find " + arguments.args.fastq_dir + " -type f -maxdepth 1 -name \'*.fastq.gz\' | parallel --jobs 8 bash " + os.path.dirname(os.path.abspath(__file__)) + "/../bash/count_fastq_reads.sh {} \'>>\' " + os.path.dirname(os.path.abspath(__file__)) + "/../../output/fastq_read_counts.txt"
	os.system(bash_command)
	#Read in outputted reads text file add column names to dataframe
	reads_df = pd.read_csv(os.path.dirname(os.path.abspath(__file__)) + "/../../output/fastq_read_counts.txt", sep="\s+", names=["Sample", "Reads"])
	#Save dataframe to file
	reads_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + "/../../output/reads_file.tsv", sep="\t", index=False)
	os.remove(os.path.dirname(os.path.abspath(__file__)) + "/../../output/fastq_read_counts.txt")
	print("Reads file has been created!")
	print("WARNING: Path names in script must be edited to represent Sample names!")


def create_binsize_file(arguments):
	"""
	This function creates the binsize file needed.
	Input(s):
	arguments are the command line arguments.
	Output(s):
	A file containing the bin size information is saved in the "output" folder.
	"""

	#Run bach script to create a tab file containing the bin size information
	bash_command = "sh " + os.path.dirname(os.path.abspath(__file__)) + "/../bash/calc_genome_size.sh " + arguments.args.fna_dir
	os.system(bash_command)
	#Read in  the outputted binsize tab file and add column names
	binsize_df = pd.read_csv(arguments.args.fna_dir + "genomesize.tab", sep="\t", names=["Bin", "Size"])
	#Save dataframe to file
	binsize_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + "/../../output/binsize_file.tsv", sep="\t", index=False)
	#Delete unneeded tab file
	os.remove(arguments.args.fna_dir + "genomesize.tab")
	#Delete files created from bash script
	for bins in glob.glob(arguments.args.fna_dir + '*.totalbp.txt'):
		os.remove(bins)
	print("Bin size file has been created!")


def create_depth_file(arguments, depth_format, depth_dir):
	"""
	This function creates the depth file needed.
	Input(s):
	arguments are the command line arguments.
	depth_df is a dataframe containing the depth file in a wide format.
	Output(s):
	A file containing the depth information is saved in the "output" folder.
	"""

	#If depth txt files are not concatenated already
	if depth_dir:
		depth_list_1 = []
		depth_list_2 = []
		depth_concat_list = []

		#Format directory
		if depth_format[-1] != '/':
			depth_format = depth_format + '/'
		#Iterate through depth files in directory
		for file in glob.glob(depth_format+'*'):
			if 'depth' in file:
				if 'txt' in file:
					depth_df = pd.read_csv(file, sep='\s+')
				elif 'tsv' in file:
					depth_df = pd.read_csv(file, sep='\t')
				elif 'csv' in file:
					depth_df = pd.read_csv(file)
				##If column exists, drop it
				if 'totalAvgDepth' in depth_df.columns:
					depth_df = depth_df.drop(columns=['totalAvgDepth'])
				#Remove all columns containing var in them
				depth_df = depth_df[depth_df.columns.drop(list(depth_df.filter(regex='var')))]
				#If each bin was mapped to all assemblies
				if len(depth_df.columns) > 3:
					depth_df = format_depth_df(depth_df)
					depth_list_1.append(depth_df)
				else:
					depth_list_2.append(depth_df)

		if len(depth_list_2) > 0:
			i = 1
			for df in depth_list_2:
				if i == 1:
					depth_df = copy.deepcopy(df)
				else:
					depth_df = depth_df.merge(df, on=['contigName', 'contigLen'], how='outer')
				i+=1
			depth_df = format_depth_df(depth_df)
			depth_concat_list.append(depth_df)
		if len(depth_list_1) > 0:
			depth_df = pd.concat(depth_list_1)
			depth_concat_list.append(depth_df)
		if len(depth_list_1) > 0 and len(depth_list_2) > 0:
			depth_df = pd.concat(depth_concat_list)
	else:
		depth_df = depth_format
		if 'totalAvgDepth' in depth_df.columns:
			depth_df = depth_df.drop(columns=['totalAvgDepth'])
		depth_df = depth_df[depth_df.columns.drop(list(depth_df.filter(regex='var')))]
		depth_df = format_depth_df(depth_df)

	#Save dataframe to file
	depth_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + "/../../output/depth_file.tsv", sep="\t", index=False)
	print("Depth file has been created! Samples containing no depth information were removed.")
	print("WARNING: If sample names in the \"Sample\" column does not match the other input files, it must be manually edited!")


def format_depth_df(depth_df):
	"""
	This function formats the depth dataframe accordingly.
	Input(s):
	depth_df is a pandas dataframe that contains all the depth information for each of the bins.
	Output(s):
	depth_df is a pandas dataframe that contains all the depth information for each of the bins.
	"""

	#Rename certain columns
	depth_df = depth_df.rename(columns=lambda x: re.sub('_S\d+','',x))
	#Pivot dataframe into a long format
	depth_df = pd.melt(depth_df, id_vars=['contigName', 'contigLen'], value_vars=depth_df.columns.tolist()[2:], var_name='Sample',value_name='Depth')
	#Rename a column
	depth_df = depth_df.rename(columns={"contigName": "Original_Contig_Name"})
	depth_df = depth_df.dropna()

	return depth_df


def create_mapping_file(arguments, bin_sample):
	"""
	This function creates the mapping file needed.
	Input(s):
	arguments are the command line arguments.
	bin_sample is a dataframe mapping each bin to its respective sample.
	Output(s):
	A file containing the mapping information is saved in the "output" folder.
	"""

	contig_lst = []
	og_contig_lst = []
	bin_lst = []

	#Get contig names
	for file in glob.glob(arguments.args.fna_dir+"*.fna"):
		f = open(file, "r")
		lines = f.readlines()
		for l in lines:
			if ">" in l:
				contig_lst.append(l[1:-1])
				#Get bin names for mapping
				bin_lst.append(os.path.basename(file)[:-4])

	#Add each contig name to each bin name in bin_sample
	mapping_df = pd.DataFrame(zip(contig_lst, bin_lst), columns=["Original_Contig_Name", "Bin"])
	mapping_df = mapping_df.merge(bin_sample, on="Bin", how="left")
	mapping_df = mapping_df.dropna()

	#Edit names in Original_Contig_Name column
	i = 0
	for sample in mapping_df['Original_Contig_Name']:
		idx = re.search("scaffold", sample)
		if idx != None:
			og_contig_lst.append(mapping_df['Sampling_Site'].iloc[i]+'_'+sample[idx.start():])
		else:
			og_contig_lst.append(np.nan)
		i+=1
	mapping_df['Original_Contig_Name'] = og_contig_lst
	mapping_df = mapping_df.dropna()
	#Save dataframe to file
	mapping_df.to_csv(os.path.dirname(os.path.abspath(__file__)) + "/../../output/mapping_file.tsv", sep="\t", index=False)
	print("Mapping file has been created!")


def main():

	arguments = Command_line_args()

	#Create output directory if not already present
	if not os.path.exists(os.path.dirname(os.path.abspath(__file__)) + "/../../output"):
		os.makedirs(os.path.dirname(os.path.abspath(__file__)) + "/../../output")

	#Add backslash to the end of path
	if arguments.args.fastq_dir[-1] != "/":
		arguments.args.fastq_dir = arguments.args.fastq_dir + "/"
	#Add backslash to the end of path
	if arguments.args.fna_dir[-1] != "/":
		arguments.args.fna_dir = arguments.args.fna_dir + "/"

	#Check if files are tsv or csv
	if arguments.args.bin_samples.endswith('.tsv'):
		bin_sample_df = pd.read_csv(arguments.args.bin_samples, sep="\t")
	elif arguments.args.bin_samples.endswith('.csv'):
		bin_sample_df = pd.read_csv(arguments.args.bin_samples)
	elif arguments.args.bin_samples.endswith('.txt'):
		bin_sample_df = pd.read_csv(arguments.args.bin_samples, sep='\s+')
	else:
		print("Bin to sample file is not in tsv or csv format!")
		quit()

	#Create depth file
	print("Beginning to create input files. This may take a while.")
	if arguments.args.depth:
		if arguments.args.depth.endswith('.tsv'):
			depth_df = pd.read_csv(arguments.args.depth, sep="\t")
		elif arguments.args.depth.endswith('.csv'):
			depth_df = pd.read_csv(arguments.args.depth)
		elif arguments.args.depth.endswith('.txt'):
			depth_df = pd.read_csv(arguments.args.depth, sep='\s+')
		else:
			print("Depth file is not in tsv, csv, or txt format!")
			quit()
		create_depth_file(arguments, depth_df, False)
	elif arguments.args.depth_dir:
		create_depth_file(arguments, arguments.args.depth_dir, True)
	else:
		print("Must have the depth or depth_dir arguments! Please rerun with the required arguments.")
		quit()

	#Create input files
	create_reads_file(arguments)
	create_binsize_file(arguments)
	create_mapping_file(arguments, bin_sample_df)
	print("Success!\nThe following files have been saved in the \"output\" directory:\n\ndepth_file.tsv\nreads_file.tsv\nbinsize_file.tsv\nmapping_file.tsv\n")
	print("WARNING: Verify that all the scaffold names in the depth and mapping files under the \"Original_Contig_Name\" column are formatted the same. If not, double check the headers in the initial depth file(s) used.\n")


if __name__ == "__main__":
	main()
