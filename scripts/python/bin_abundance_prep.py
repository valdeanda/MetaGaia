#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     bin_abundance_prep.py
# Purpose:  create the input files needed to calculate bin abundacnces.
#
# @uthors:      vda - valdeanda@utexas.edu, acph  dragopoot@gmail.com
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

def create_reads_df(arguments):
	"""

	"""

	bash_command = "find " + arguments.args.fastq_dir + " -type f -name \'*.fastq.gz\' | parallel --jobs 8 bash ../bash/count_fastq_reads.sh {} \'>>\' ../../output/fastq_read_counts.txt"
	os.system(bash_command)
	reads_df = pd.read_csv("../../output/fastq_read_counts.txt", sep="\s+")
	reads_df.columns = ["Sample", "Reads"]
	reads_df.to_csv("../../output/reads_file.tsv", sep="\t", index=False)
	print("Reads file has been created!")
	print("WARNING: Path names in script must be edited to represent sample names!")

def create_binsize_df(arguments):
	"""

	"""

	bash_command = "sh ../bash/calc_genome_size.sh " + arguments.args.fna_dir
	os.system(bash_command)
	binsize_df = pd.read_csv(arguments.args.fna_dir + "genomesize.tab", sep="\t", names=["Bin", "Size"])
	binsize_df.to_csv("../../output/binsize_file.tsv", sep="\t", index=False)
	os.remove(arguments.args.fna_dir + "genomesize.tab")
	print("Bin size file has been created!")

def create_depth_df(arguments, depth):
	"""

	"""

	print("WARNING: Please make sure the columns of the depth file represent the sample names and not bam file names!")

	depth_df = pd.read_csv(depth, sep="\t")
	depth_df = depth_df.unstack().reset_index()
	depth_df.rename({"level_0": "Sample", 0: "Sample_Depth", "contigName": "Original_Contig_Name"})
	depth_df.to_csv("../../output/depth_file.tsv", sep="\t", index=False)
	print("Depth file has been created!")

def create_mapping_df(arguments, bin_sample):
	"""

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
	if ".tsv" in bin_samples:
		bin_sample_df = pd.read_csv(bin_samples, sep="\t")
	elif ".csv" in bin_samples:
		bin_sample_df = pd.read_csv(bin_samples)
	else:
		print("Bin to sample file is not in tsv or csv format!")
		quit()
	if ".tsv" in depth:
		depth_df = pd.read_csv(depth, sep="\t")
	elif ".csv" in depth:
		depth_df = pd.read_csv(depth)
	else:
		print("Depth file is not in tsv or csv format!")
		quit()

	#Create input files
	print("Beginning to create input files. This may take a while.")
	create_reads_df(arguments)
	create_binsize_df(arguments)
	create_depth_df(arguments, depth_df)
	create_mapping_df(arguments, bin_sample_df)
	print("All input files have been successfully created and can be found in the \"output\" directory!")


if __name__ == "__main__":
	main()

















