#!/usr/bin/env python

# This script represents a web crawler for www.uniprot.org.
#
# This crawler is able to fetch the content
# for a number of given protein ids (or p_id) and parse the kind and number of
# mutations out of it.
#
# The fetched mutations will then be processed and shown to
# the user in form of:
# - a list of mutations printed to the console
# - a list of mutations saved within a csv file
# - a list of the fetched proteins with their corresponding name
#   within a csv file (meta information)
#
# For every call of this script, it will create a directory where
# the result of this crawl is saved to.
# The can be manually done by prepending a header line within the passed input file.
# If no header line is given, a unique UUID will be generated.
#
# [OPTIONAL] Futhermore, this crawler is also able present the
# the fetched mutations in a with matplotlib generated colormap. This colormap
# This features resembles the contribution chart on Github (https://goo.gl/images/Hb8YoE).
# Please note, that in order to use this feature, the necessary packages must be
# installed beforehand. (see requirements.txt for the individual packages).
# The best way to install the packages would be with pip:
# $ pip install -r requirements.txt
#
# This script was developed with Python 3.6.2

# General imports
import re
import itertools
import math
import sys
import os
import uuid
import urllib.request
import socket
from functools import reduce

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

# Conducts the crawling to uniprot.org and parses the transmitted content
# into a mutation dictionary and a protein list
def crawl(app):
	print("[i] reading uniprot ids from {0}".format(app["input_file"]))

	with open(app["input_file"]) as filehandler:
		# a list of protein names
		proteins = []

		# contains all mutatations as a nested dictionary
		mutations = {}

		# iterate over ids
		for p_id in filehandler:

			# sanitize input
			p_id = p_id.strip()

			# skip header line
			if p_id.startswith(">"):
				continue

			target = os.path.join(UNIPROT_URL, p_id)
			print("[i] crawling {0}".format(target))
			try:
				# send HTTP GET request
				request = urllib.request.Request(url=target, method="GET")
				response = urllib.request.urlopen(request, timeout=app["timeout"])
				if response.status != 200:
					print("[x] failed with {0}: {1}".format(response.status, response.reason))
					continue

				# decode content
				content = response.read().decode("utf-8")

			except (urllib.error.HTTPError, urllib.error.URLError):
				print("[x] could not connect to {0}".format(target))
				continue

			except socket.timeout:
				print("[x] could not connect to {0}: host does not respond".format(target))
				continue

			# retrive dictionary of variants
			protein, mutations = retrieve_variants(content, mutations)

			if protein is "":
				print("[x] no mutations found for {0}, skipped...".format(target))
				continue

			proteins.append({
				"id": p_id,
				"name": protein
			})

	return proteins, mutations

# Parses variants out of the given content into a mutations dictionary
def retrieve_variants(content, mutations):
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
def prepare_statistic_data(mutations):

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
def draw_statistic(mutations, app):

	# Graph-specific imports
	import matplotlib.pylab as pl
	import numpy as np
	import pandas as pd

	# For serializing the plot
	import dill

	# if no mutations, return
	if not mutations:
		print("[w] no mutations given. cannot create graph")
		return

	columns, rows, values, raw = prepare_statistic_data(mutations)

	# create a R-like dataframe
	df = pd.DataFrame({"fromAA": columns, "toAA": rows, "weight": values})

	# reshape the data
	dataframe = df.pivot(index="fromAA", columns="toAA", values="weight")

	# create a meshgrid for the column and rows
	Rows, Columns = np.mgrid[:dataframe.shape[0]+1, :dataframe.shape[1]+1]

	# create subplots
	fig, ax = pl.subplots(figsize=(9, 6))
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
	pl.title("Mutations in {0}".format(app["name"]), y=1.06, fontsize=10)
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
	figurepic_path = os.path.join(app["path"], "variants_mesh.png")
	fig.savefig(figurepic_path)
	print("[i] figure picture saved to {0}".format(figurepic_path))

	# save the figure as binary file to reopen it later
	figurebin_path = os.path.join(app["path"], "variants_mesh.matpltlib")
	with open(figurebin_path, "wb") as filehandler:
		dill.dump(ax, filehandler)
		print("[i] figure binary saved to {0}".format(figurebin_path))

# Print all given mutations of all proteins to the console
def print_to_console(mutations, app):
	print("Mutations in {0}".format(app["name"]))
	print("="*(13 + len(app["name"])))

	for key, variants in sorted(mutations.items()):
		print("[Key]: {0}".format(key))
		for variant, number in sorted(variants.items()):
			print(" {0} {1} {2} ({3})".format(" "*5, ">", variant, number))

		print()

# Writes all given mutations of all proteins to a .csv file.
def write_result_to_csv(mutations, app):
	path = os.path.join(app["path"], "variants.csv")
	with open(path, "w") as filehandler:
		filehandler.write("Mutations in {0}:\n".format(app["name"]))
		filehandler.write("from AA,to AA,number of variants\n")

		for key, variants in sorted(mutations.items()):
			for variant, number in sorted(variants.items()):
				filehandler.write("{0},{1},{2}\n".format(key, variant, number))

	print("[i] csv saved to {0}".format(path))

# Writes all used proteins with id and name to a .csv file.
# This represents meta information
def write_meta_to_csv(proteins, app):
	path = os.path.join(app["path"], "meta.csv")
	with open(path, "w") as filehandler:
		filehandler.write("Crawled proteins in {0}:\n".format(app["name"]))
		filehandler.write("id,name\n")

		for protein in proteins:
			filehandler.write("{0},{1}\n".format(protein["id"], protein["name"]))

	print("[i] meta data saved to {0}".format(path))

# Parses optional command line arguments
# and set the fields for the application
def parse_cmd_args(app):
	# program parameter passed by cmd
	args = sys.argv

	app["input_file"] = args[1] if len(args) > 1 else app["input_file"]
	app["output_dir"] = args[2] if len(args) > 2 else app["output_dir"]
	app["graph"] = True if len(args) > 3 and args[3] == "-graph" else app["graph"]
	app["timeout"] = int(args[4]) if len(args) > 4 else app["timeout"]

# Validates the previously parsed application arguments
# and does some preparation
def validate_cmd_args(app):
	# check whether input file exists
	if not os.path.exists(app["input_file"]):
		return "input file {0} does not exist. Abort..".format(app["input_file"]), True

	# check whether input file is empty
	if os.stat(app["input_file"]).st_size == 0:
		return "input file {0} is empty. Abort..".format(app["input_file"]), True

	# create output directory
	os.makedirs(app["output_dir"], exist_ok=True)

	# retrieve name of the application
	# search for a header line in input file
	# if no one is given, define unique uuid as name
	app["name"] = str(uuid.uuid4())

	# check whether a header line is existing
	with open(app["input_file"]) as filehandler:
		header = filehandler.readlines()[0]

		if header.startswith(">"):
			# get and sanitize the name
			sanitized = header[1:]
			sanitized = sanitized.strip()
			sanitized = sanitized.replace(" ", "_")

			# only take the sanitized head line if it's not empty
			if sanitized is not "":
				app["name"] = sanitized

	# set joined result path where the fetched result will be saved to
	app["path"] = os.path.join(app["output_dir"], app["name"])

	# create result directory
	os.makedirs(app["path"], exist_ok=True)

	return "", False

# The main method in which the program runs
def main():
	# define application, default program parameter
	app = {
		"input_file": "uniprot_ids",
		"output_dir": "output/",
		"graph": False,
		"timeout": 15,
	}

	# parse command line arguments
	parse_cmd_args(app)

	# validate command line arguments
	message, err = validate_cmd_args(app)
	if err:
		print("[x] {0}".format(message))
		sys.exit(1)

	# start crawling
	print("[i] starting uniprot crawler")
	proteins, mutations = crawl(app)

	# print the data to console
	print_to_console(mutations, app)

	# save meta data into a .csv file
	write_meta_to_csv(proteins, app)

	# save result data (mutations) into a .csv file
	write_result_to_csv(mutations, app)

	# draw the statistic
	if app["graph"]:
		draw_statistic(mutations, app)

	# stop crawling
	print("[i] stopping uniprot crawler")

if __name__ == '__main__':
	main()