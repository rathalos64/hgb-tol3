#!/usr/bin/env python

# General imports
import re
import itertools
import math
import sys
import os
import uuid
from functools import reduce
import requests

# The Uniprot link
UNIPROT_URL = "http://www.uniprot.org/uniprot/"

# Needed regular expressions for parsing the HTTP content
# u"\u2192" == "→"
REGEX_MUTATIONS = "{0} {1} {2}".format(
	r"<a href=\"(?:.+?)\">(?:(\w)",
	u"\u2192",																																																	
	r"(\w))</a>"
)
REGEX_NAME = "<h1 property=\"schema\:name\">(.+?)<\/h1>"

# parses variants out of the given content into a mutations dictionary
def retrieveVariants(content, mutations):
	# parse mutations
	for variants in re.compile(REGEX_MUTATIONS).findall(content):
		fromA = variants[0]
		toA = variants[1]

		if fromA not in mutations.keys():
			mutations[fromA] = {}

		if toA not in mutations[fromA].keys():
			mutations[fromA][toA] = 0

		mutations[fromA][toA] = mutations[fromA][toA] + 1

	# parse name
	name = re.compile(REGEX_NAME).search(content).group(1)
	return name, mutations

# Prepare the mutations to align with a pcolor diagram in matplotlib.
def prepereDataForStatistic(mutations):

	####################### [COLUMNS] #######################

	columns = []
	for key in mutations.keys():
		for i in range(0, len(mutations.keys())):
			columns.append(key)
	columns = sorted(columns)

	######################### [ROWS] #########################

	rows = sorted(list(set(
		reduce(lambda acc, curr: acc + curr,
			[tuple(row) for row in map(lambda x: x.keys(), mutations.values())],
			()
		)
	)))
	n = len(rows)
	tmp = []
	for i in range(0, n):
		tmp.extend(rows)
	rows = tmp

	######################## [VALUES] #######################

	values = []
	for key in sorted(set(columns)):
		entry = []
		for value in sorted(set(rows)):
			x = mutations[key][value] if value in mutations[key].keys() else 0
			entry.append(x)
		values.append(entry)

	raw = values
	values = list(itertools.chain.from_iterable(values))
	
	return columns, rows, values, raw

# Draw the statistic with matplotlib by using a pcolor / pcolormesh representation.
def drawStatistic(proteins, mutations, base_path):

	# Graph-specific imports
	import matplotlib.pylab as pl
	import numpy as np
	import pandas as pd

	# For serializing the plot
	import dill

	columns, rows, values, raw = prepereDataForStatistic(mutations)

	# create a R-like dataframe
	df = pd.DataFrame({"fromAA": columns, "toAA": rows, "weight": values})

	# reshape the data
	dataframe = df.pivot(index="fromAA", columns="toAA", values="weight")

	# create a meshgrid for the column and rows
	Rows, Columns = np.mgrid[:dataframe.shape[0]+1, :dataframe.shape[1]+1]

	# create subplots
	fig, ax = pl.subplots(figsize=(7, 5))
	ax.set_aspect("equal")
	ax.xaxis.set_ticks_position('top')

	# draw the pcolor graph
	pl.pcolor(Columns, Rows, dataframe.values, cmap="Greens", edgecolor="w", vmin=0, vmax=max(values))

	# set options
	rowsX = sorted(set(rows))
	columnsY = sorted(set(columns))
	pl.xticks(range(len(rowsX)), rowsX, fontsize=8, rotation=45)
	pl.yticks(range(len(columnsY)), columnsY, fontsize=8, rotation=45)
	pl.ylim(max(pl.ylim()), min(pl.ylim()))

	# set labels
	pl.ylabel("from amino acid")
	pl.xlabel("to amino acid")
	pl.title("Mutations in {0}".format(proteins), y=1.06, fontsize=10)
	pl.tight_layout(pad=1.5)

	# set axis options
	ax.tick_params(axis='both', which='major', pad=2)
	ax.patch.set_facecolor('grey')

	ax.xaxis.set(ticks=np.arange(0.5, len(rowsX)), ticklabels=rowsX)
	ax.yaxis.set(ticks=np.arange(0.5, len(columnsY)), ticklabels=columnsY)

	# the right legend color bar
	pl.colorbar()

	# define a function for setting the left lower corner tooltip values 
	def format_coord(x, y):
		row = math.floor(y)
		col = math.floor(x)
		if row < 0 or row > len(raw) - 1:
			return ""

		if col < 0 or col > len(raw[row]) - 1:
			return ""

		return "{0} {1} {2} | variants = {3}".format(sorted(set(columns))[row], ">", sorted(set(rows))[col], raw[row][col])

	ax.format_coord = format_coord

	# save the figure as a picture
	figurepic_path = os.path.join(base_path, "variants_mesh.png")
	fig.savefig(figurepic_path)
	fig.canvas.set_window_title('Mutations (natural variants) in proteins')
	print("[i] figure picture saved to '{0}'".format(figurepic_path))

	# save the figure as binary file to reopen it later
	figurebin_path = os.path.join(base_path, "variants_mesh.matpltlib")
	with open(figurebin_path, "wb") as filehandler:
		dill.dump(ax, filehandler)
		print("[i] figure binary saved to '{0}'".format(figurebin_path))

# Print all given mutations of all proteins to the console
def printToConsole(proteins, mutations):
	print("Mutations in {0}".format(proteins))
	print("="*25)

	for key, variants in sorted(mutations.items()):
		print("[Key]: {0}".format(key))
		for variant, number in sorted(variants.items()):
			print(" {0} {1} {2} ({3})".format(" "*5, "-->", variant, number))

		print()

# Write all given mutations of all proteins to a .csv file.
# The path of the CSV file is given.
def printToCSV(proteins, mutations, base_path):
	path = os.path.join(base_path, "variants.csv")
	with open(path, "w") as filehandler:
		filehandler.write("Mutations in {0}:\n".format(proteins))
		filehandler.write("from AA,to AA,number of variants\n")

		for key, variants in sorted(mutations.items()):
			for variant, number in sorted(variants.items()):
				filehandler.write("{0},{1},{2}\n".format(key, variant, number))

	print("[i] csv saved to '{0}'".format(path))

def prepareResult(input_file, output_path):
	# the default name for the result directory
	# is an unique uuid
	name = str(uuid.uuid4())

	# check whether a header line is existing
	with open(input_file) as filehandler:
		header = filehandler.readlines()[0]

		if header.startswith(">"):
			# get and sanitize the name
			name = header[1:]
			name = name.strip()
			name = name.replace(" ", "_")

	# final path for saving the result
	path = os.path.join(output_path, name)

	# create the necessary directory
	os.makedirs(path, exist_ok=True)

	return path

# The main method in which the program runs
def main():
	# program parameter
	input_file = "proteins"
	output_path = "output/"
	with_graph = False

	# check whether file is empty or does not exist
	if os.stat(input_file).st_size == 0:
		print("[x] path {0} is either empty or does not exist. Abort..".format(input_file))
		sys.exit(1)

	# create output directory
	os.makedirs(output_path, exist_ok=True)

	# start crawling
	print("[i] starting uniprot crawler")
	print("[i] reading uniprot ids from '{0}'".format(input_file))

	with open(input_file) as filehandler:
		# a list of protein names
		proteins = []

		# contains all mutatations as a nested dictionary
		mutations = {}

		# iterate over ids
		for p_id in filehandler:

			# sanitize input
			p_id = p_id.strip()

			# skip proteins header
			if p_id.startswith(">"):
				continue

			# send HTTP GET request
			print("[i] crawling for '{0}'".format(p_id))

			resp = requests.get(os.path.join(UNIPROT_URL, p_id))
			if resp.status_code != 200:
				print("[x] failed with {0}: {1}".format(resp.status_code, resp.reason))
				print("[i] skip")
				continue

			content = resp.content.decode("utf-8")

			# retrive dictionary of variants
			protein, mutations = retrieveVariants(content, mutations)

			if protein is "":
				print("[x] no mutations found for '{0}', skipped...".format(p_id))
				continue

			proteins.append(protein)

	# prepare the result of the request
	path = prepareResult(input_file, output_path)

	# print the data to console
	printToConsole(proteins, mutations)

	# save the data into a .csv file
	printToCSV(proteins, mutations, path)

	# draw the statistic
	if with_graph:
		drawStatistic(proteins, mutations, path)

if __name__ == '__main__':
	main()