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
import functools
import glob
import math
import numpy as np
import os
import pandas as pd
import random
import subprocess

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
		#self.parser.add_argument('bin_samples', type=str, help='Input file path containing each bin mapped to each sample (with extension).')
		#self.parser.add_argument('depth', type=str, help='Input file path containing depth information (with extension).')
		self.parser.add_argument('fasta_dir', type=str, help='Input directory path containing all the fasta files.')
		self.parser.add_argument('fna_dir', type=str, help='Input directory path containing all the fna files.')
		self.args = self.parser.parse_args()

def create_reads_df(arguments):
	"""

	"""

	subprocess.call(['bash', '../bash/count_fastq_reads.sh', arguments.args.fasta_dir])
	for file in os.listdir(arguments.args.fasta_dir):
		if file.startswith("fastq_read_counts"):
			reads_df = pd.read_csv(arguments.args.fasta_dir + str(file), sep="\t")
	reads_df.columns = ["Sample", "Reads"]
	reads_df.to_csv("../../output/reads_file.tsv", index=False)
	os.remove(arguments.args.fasta_dir + str(file))
	print("Reads file has been created!")
	print("WARNING: Path names in script must be edited to represent sample names!")

def create_binsize_df(arguments):
	"""

	"""

	subprocess.call(['bash', '../bash/calc_genome_size.sh', arguments.args.fna_dir])
	binsize_df = pd.read_csv(arguments.args.fna_fir + "genomesize.tab", sep="\t")
	binsize_df.columns = ["Bin", "Size"]
	binsize_df.to_csv("../../output/binsize_file.tsv", index=False)
	os.remove(arguments.args.fasta_dir + "genomesize.tab")
	print("Bin size file has been created!")

# def create_depth_df(arguments):
# 	"""

# 	"""

# 	depth_df = pd.read_csv(depth, sep="\t")
# 	depth_df = pd.melt(depth_df, id_vars=)


def main():

	arguments = Command_line_args()

	#Add backslash to the end of path
	if arguments.args.fasta_dir[-1] != "/":
		arguments.args.fasta_dir = arguments.args.fasta_dir + "/"
	#Add backslash to the end of path
	if arguments.args.fna_dir[-1] != "/":
		arguments.args.fna_dir = arguments.args.fna_dir + "/"

	#Create input files
	print("Beginning to create input files.")
	reads_df = create_reads_df(arguments)
	binsize_df = create_binsize_df(arguments)
	# depth_df = create_depth_df(arguments)
	# mapping_df = create_mapping_df(arguments)
	print("All input files have been successfully created and can be found in the \"output\" directory!")


if __name__ == "__main__":
	main()

















