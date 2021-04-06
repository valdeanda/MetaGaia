#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# ------------------------------
# Name:     metagaia_viz.py
# Purpose:  use the output file after calculating bin abundacnces to visualize the data.
#
# @uthors:      vda - valdeanda@utexas.edu, acph  dragopoot@gmail.com
#
# Created:     2021
# Licence:     GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
# ------------------------------
import argparse
import functools
import math
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import re
import seaborn as sns

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
		self.parser.add_argument('bin_abundance', type=str, help='Input file path outputted from MetaGaia with bin abundances with extension (str).')
		self.parser.add_argument('taxonomy_info', type=str, help='Input file path containing taxonomy information (with extension).')
		self.parser.add_argument('percent', type=float, nargs='?', default=10, help='Percent of highest sample in each bin [10] (float).')
		self.parser.add_argument('width', type=int, nargs='?', default=4, help='Width of outputted clustermap figure [4] (int).')
		self.parser.add_argument('height', type=int, nargs='?', default=5, help='Height of outputted clustermap figure [5] (int).')
		self.parser.add_argument('dpi', type=int, nargs='?', default=300, help='Resolution for output figure file [300] (int).')
		self.parser.add_argument('out_fig', type=str, nargs='?', default="test.png", help='Stores the figure in the specified file path and format [test.png] (str).')
		self.parser.add_argument('taxa_color', type=str, default="", nargs='?', help='Input file path containing the color code for each taxa with extension [""] (str).')
		self.args = self.parser.parse_args()


def format_dataframe(arguments, bin_abundances, taxonomy_info):
	"""
	Formats the bin abundances dataframe to a layout useable for cluster mapping.
	Input(s):
	arguments is a class containing all the command line arguments.
	bin_abundances is a pandas dataframe containing all the bin abundance information.
	taxonomy_info is a pandas dataframe containing all the taxaonomy information for each bin.
	Output(s):
	A formatted bin abundances dataframe.
	"""

	#Set Sample column to the site number
	for sample in bin_abundances['Sample'].unique():
		digit_index = re.search("_[0-9]", sample)
		bin_abundances['Sample'] = bin_abundances.Sample.str.replace(sample, sample[digit_index.start()+1:])

    #Drop unneeded columns
	bin_abundances = bin_abundances.drop(columns=['RelativeAbundance'])
	#Log abundances
	bin_abundances['RelativeAbundanceReadable'] = np.log10(bin_abundances['RelativeAbundanceReadable'])
    #Pivot table wider so that each bin is mapped to its respective site
	bin_abundances = bin_abundances.pivot(index='Bin', columns='Sample', values='RelativeAbundanceReadable')
    #Replace NaN with 0 and set the index to the Bin
	bin_abundances = bin_abundances.replace(np.nan, 0)

	if taxonomy_info.columns.tolist() != ['Bin', 'Taxa']:
		#Rename headers
		taxonomy_info.columns = ['Bin', 'Taxa']
	#Left join bin abundance file with taxanomy information
	bin_abundances = bin_abundances.merge(taxonomy_info, on='Bin', how='left').set_index('Bin')
	#Order dataframe by Taxa
	bin_abundances = bin_abundances.sort_values('Taxa')
	#Sort columns
	bin_abundances = bin_abundances.reindex(sorted(bin_abundances.columns), axis=1)
	#Replace underscore with space in taxa names
	bin_abundances = bin_abundances.replace('_', ' ', regex=True)

	return bin_abundances


def get_color_palette(bin_abundances):
	"""
	Obtains the number of colors needed to label each taxa and converts RGB to hex codes.
	Input(s):
	bin_abundances is a pandas dataframe with all the information joined together.
	Output(s):
	A list containing colors for heatmap.
	"""

	#Fix how color scheme is selected
	final_colors_list = []

	#Get hex color codes
	for col in range(len(bin_abundances['Taxa'].unique())):
		r = lambda: random.randint(0,255)
		color = '#%02X%02X%02X' % (r(),r(),r())
		while color in final_colors_list:
			r = lambda: random.randint(0,255)
			color = '#%02X%02X%02X' % (r(),r(),r())
		final_colors_list.append(color)

	return final_colors_list


def filter_top_sites(arguments, bin_abundances):
	"""
	Creates a dataframe with the top taxa in each bin.
	Input(s):
	arguments is a class containing all the command line arguments.
	bin_abundances is a pandas dataframe with all the information joined together.
	Output(s):
	merged_dfs is a new pandas dataframe containing the top samples in each column.
	"""

	#List to keep dataframes together
	df_list = []
	#Get list of columns
	col_list = bin_abundances.drop(columns = 'Taxa').columns.tolist()
	#Number of rows to keep for each bin
	num_filter = math.floor((arguments.args.percent / 100) * len(bin_abundances))

	#Create new dataframe with top values for each column
	for col in col_list:
		#List of columns to drop after getting rows to keep
		drop_list = [bins for bins in col_list if bins != col]
		#Get needed rows
		bin_abundances_cp = bin_abundances.nlargest(n = num_filter, columns = col)
		#Drop all columns not of interest
		bin_abundances_cp = bin_abundances_cp.drop(columns = drop_list)
		df_list.append(bin_abundances_cp)

	#Combine dataframes into one
	merged_dfs = pd.concat(df_list)
	merged_dfs = merged_dfs.replace(np.nan, 0)
	merged_dfs = merged_dfs.sort_values('Taxa')
	merged_dfs = merged_dfs.drop_duplicates()
	merged_dfs = merged_dfs.reindex(sorted(merged_dfs.columns), axis=1)

	return merged_dfs


def create_clustermap(arguments, bin_abundances, palette, row_colors, heatmp = False):
	"""
	Creates clustermap with taxonomic information.
	Input(s):
	bin_abundances is a pandas dataframe with all the information joined together.
	palette is a dictionary containing the taxa-to-color information.
	row_colors is a pandas object that maps the taxa to the color.
	Output(s):
	A clustermap of all the bins and taxa.
	"""

	#Get the greatest dimension value
	max_size = max(arguments.args.width, arguments.args.height)

	#Set font size
	sns.set(font_scale=max_size/10)

	#Handles for taxa color legend
	handles = [Patch(facecolor=palette[name]) for name in palette]

	#Remove y axis labels if there are too many
	if len(bin_abundances) > 50:
		y_tick_label = False
	else:
		y_tick_label = True
	#Create clustermap with certain specifications
	if heatmp:
		cluster = sns.clustermap(bin_abundances.drop(columns="Taxa"), col_cluster=False, row_cluster=False, yticklabels=y_tick_label, cmap="Reds", row_colors=row_colors, method='ward', figsize=(arguments.args.width, arguments.args.height), cbar_pos=(0.05, .25, .03, .4))
		#Create taxa-color legend
		plt.legend(handles, palette, bbox_to_anchor=(1.2, 0.82), bbox_transform=plt.gcf().transFigure, loc='upper right', prop={'size': max_size/2}, facecolor="white", edgecolor="white")
	else:
		cluster = sns.clustermap(bin_abundances.drop(columns="Taxa"), col_cluster=True, row_cluster=True, yticklabels=y_tick_label, cmap="Reds", row_colors=row_colors, method='ward', figsize=(arguments.args.width, arguments.args.height))
		#Create taxa-color legend
		plt.legend(handles, palette, bbox_to_anchor=(1.2, 1), bbox_transform=plt.gcf().transFigure, loc='upper right', prop={'size': max_size/2}, facecolor="white", edgecolor="white")

	#Save or show plot
	if arguments.args.out_fig:
		if heatmp:
			plt.savefig("../../output/heatmap_top_" + arguments.args.out_fig, dpi=arguments.args.dpi, bbox_inches="tight")
		else:
			plt.savefig("../../output/clustermap_" + arguments.args.out_fig, dpi=arguments.args.dpi, bbox_inches="tight")
	else:
		plt.show()


def main():

	arguments = Command_line_args()

	#Read in bin abundance file
	if "tsv" in arguments.args.bin_abundance:
		bin_abundance_df = pd.read_csv("../../data/" + arguments.args.bin_abundance, sep = "\t")
	else:
		bin_abundance_df = pd.read_csv("../../data/" + arguments.args.bin_abundance)

	#Read in taxanomy file
	if "tsv" in arguments.args.taxonomy_info:
		taxonomy_df = pd.read_csv("../../data/" + arguments.args.taxonomy_info, sep = "\t", header = None)
	else:
		taxonomy_df = pd.read_csv("../../data/" + arguments.args.taxonomy_info, header = None)
	
	#Format dataframe to desired layout
	bin_abundance_df = format_dataframe(arguments, bin_abundance_df, taxonomy_df)

    #Get colors for each taxa to use in heatmap
	if len(arguments.args.taxa_color) == 0:
		final_colors = get_color_palette(bin_abundance_df)
	else:
		if "tsv" in arguments.args.taxa_color:
			colors_df = pd.read_csv("../../data/" + arguments.args.taxa_color, sep = "\t", header = None)
		else:
			colors_df = pd.read_csv("../../data/" + arguments.args.taxa_color, header = None)
		final_colors = colors_df[0].tolist()

   	#Set color palette
	my_palette = dict(zip(bin_abundance_df.Taxa.unique(), final_colors))
	row_colors = bin_abundance_df.Taxa.map(my_palette)

	#Outputs final clustermap
	create_clustermap(arguments, bin_abundance_df, my_palette, row_colors, False)

	#Get dataframe with highest abundances in each column
	merged_df = filter_top_sites(arguments, bin_abundance_df)
	merged_df.to_csv('../../output/top_sample_abundances.tsv', sep='\t')

	#Set color palette
	my_palette = dict(zip(merged_df.Taxa.unique(), final_colors))
	row_colors = merged_df.Taxa.map(my_palette)

	#Outputs final heatmap
	create_clustermap(arguments, merged_df, my_palette, row_colors, True)

	print("Files in \'output\' folder.\nEnjoy!")

if __name__ == "__main__":
	main()

