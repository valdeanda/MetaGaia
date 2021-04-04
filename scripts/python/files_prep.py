#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     metagaia_prep.py
# Purpose:  create the input files needed to calculate bin abundacnces.
#
# @uthors:      sbs - sahil2699@gmail.com, acph  dragopoot@gmail.com
#
# Created:     2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
import argparse
import glob
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
		self.parser.add_argument('bin_samples', type=str, help='Input file path containing each bin mapped to each sample with extension. Only tsv or csv files allowed. (str)')
		self.parser.add_argument('depth', type=str, help='Input file path containing depth information (with extension). Only tsv or csv files allowed. (str)')
		self.parser.add_argument('fastq_dir', type=str, help='Input directory path containing all the fasta files. (str)')
		self.parser.add_argument('fna_dir', type=str, help='Input directory path containing all the fna files. (str)')
		self.args = self.parser.parse_args()


def create_reads_file(arguments):
	"""
	This function creates the reads file needed.
	Input(s):
	arguments are the command line arguments.
	Output(s):
	A file containing the reads information is saved in the "output" folder.
	"""

	#Run bash script to create a text file containing the reads information
	bash_command = "find " + arguments.args.fastq_dir + " -type f -name \'*.fastq.gz\' | parallel --jobs 8 bash ../bash/count_fastq_reads.sh {} \'>>\' ../../output/fastq_read_counts.txt"
	os.system(bash_command)
	#Read in outputted reads text file
	reads_df = pd.read_csv("../../output/fastq_read_counts.txt", sep="\s+")
	#Add column names to dataframe
	reads_df.columns = ["Sample", "Reads"]
	#Save dataframe to file
	reads_df.to_csv("../../output/reads_file.tsv", sep="\t", index=False)
	print("Reads file has been created!")
	print("WARNING: Path names in script must be edited to represent sample names!")


def create_binsize_file(arguments):
	"""
	This function creates the binsize file needed.
	Input(s):
	arguments are the command line arguments.
	Output(s):
	A file containing the bin size information is saved in the "output" folder.
	"""

	#Run bach script to create a tab file containing the bin size information
	bash_command = "sh ../bash/calc_genome_size.sh " + arguments.args.fna_dir
	os.system(bash_command)
	#Read in  the outputted binsize tab file and add column names
	binsize_df = pd.read_csv(arguments.args.fna_dir + "genomesize.tab", sep="\t", names=["Bin", "Size"])
	#Save dataframe to file
	binsize_df.to_csv("../../output/binsize_file.tsv", sep="\t", index=False)
	#Delete unneeded tab file
	os.remove(arguments.args.fna_dir + "genomesize.tab")
	print("Bin size file has been created!")


def create_depth_file(arguments, depth_df):
	"""
	This function creates the depth file needed.
	Input(s):
	arguments are the command line arguments.
	depth_df is a dataframe containing the depth file in a wide format.
	Output(s):
	A file containing the depth information is saved in the "output" folder.
	"""

	#Drop all columns that contain var in the names (as they are not needed)
	depth_df = depth_df[depth_df.columns.drop(list(depth_df.filter(regex='var')))]
	#Rename certain coulumns
	depth_df = depth_df.rename(columns=lambda x: re.sub('_S\d+','',x))
	#Pivot dataframe into a long format
	depth_df = pd.melt(depth_df, id_vars=['contigName', 'contigLen', 'totalAvgDepth'], value_vars=depth_df.columns.tolist()[3:], var_name='Sample_Depth',value_name='Depth')
	#Rename a column
	depth_df = depth_df.rename({"contigName": "Original_Contig_Name"})
	#Save dataframe to file
	depth_df.to_csv("../../output/depth_file.tsv", sep="\t", index=False)
	print("Depth file has been created!")
	print("WARNING: If sample names in the \"Sample_Depth\" column does not match the other input files, it must be manually edited!")


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
	bin_lst = []

	#Get contig names
	for file in glob.glob(arguments.args.fna_dir+"*.fna"):
		f = open(file, "r")
		lines = f.readlines()
		for l in lines:
			if ">" in l:
				contig_lst.append(l[1:-1])

	#Get bin names for mapping
	for c in range(len(contig_lst)):
		idx = re.search("_scaffold", contig_lst[c])
		bin_lst.append(contig_lst[c][:idx.start()])

	#Add each contig name to each bin name in bin_sample
	mapping_df = pd.DataFrame(zip(contig_lst, bin_lst), columns=["Original_Contig_Name", "Bin"])
	mapping_df = bin_sample.merge(mapping_df, on="Bin", how="left")
	#Save dataframe to file
	mapping_df.to_csv("../../output/mapping_file.tsv", sep="\t", index=False)
	print("Mapping file has been created!")


def main():

	arguments = Command_line_args()

	#Add backslash to the end of path
	if arguments.args.fastq_dir[-1] != "/":
		arguments.args.fastq_dir = arguments.args.fastq_dir + "/"
	#Add backslash to the end of path
	if arguments.args.fna_dir[-1] != "/":
		arguments.args.fna_dir = arguments.args.fna_dir + "/"

	#Check if files are tsv or csv
	if ".tsv" in arguments.args.bin_samples:
		bin_sample_df = pd.read_csv(arguments.args.bin_samples, sep="\t")
	elif ".csv" in arguments.args.bin_samples:
		bin_sample_df = pd.read_csv(arguments.args.bin_samples)
	else:
		print("Bin to sample file is not in tsv or csv format!")
		quit()
	if ".tsv" in arguments.args.depth:
		depth_df = pd.read_csv(arguments.args.depth, sep="\t")
	elif ".csv" in arguments.args.depth:
		depth_df = pd.read_csv(arguments.args.depth)
	else:
		print("Depth file is not in tsv or csv format!")
		quit()

	#Create input files
	print("Beginning to create input files. This may take a while.")
	create_reads_file(arguments)
	create_binsize_file(arguments)
	create_depth_file(arguments, depth_df)
	create_mapping_file(arguments, bin_sample_df)
	print("All input files have been successfully created and can be found in the \"output\" directory!")


if __name__ == "__main__":
	main()

