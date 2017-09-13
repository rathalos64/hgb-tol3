# General imports
import sys
import uuid
import os
import re
import socket
import urllib.request
import collections
from functools import reduce

# The unigene URL in NCBI
UNIGENE_URL = "https://www.ncbi.nlm.nih.gov/UniGene/clust.cgi"

# Needed regular expressions for parsing HTML
REGEX_ENTRYTITLE = r'<div class="EntryTitle">\n*(.+?)\n*</div>'
REGEX_ENTRYIDENT = r'<td class="EntryIdent">(.+?)</td>'
REGEX_TISSUETABLE = r'<p(?:.+?)><b>EST(?:.+?)<\/b><\/p>\n<TABLE(.+?)<\/TABLE>'
REGEX_TISSUE = r'''
				<TR(?:.+?)>\n(?:.+?) 
					<TD(?:.+?)<\/TD>\n(?:.+?) 
					<TD(?:.+?)<\/TD>\n(?:.+?) 
					<TD(?:.+?)<\/TD>\n(?:.+?)
					<TD(?:.+?)>(.+?)<\/TD>\n(?:.+?)
					<TD(?:.+?)<\/TD>\n(?:.+?)
					<TD(?:.+?)<\/TD>\n
				<\/TR>'''

# Conducts the crawling to the unigene database and parses the transmitted content
# into a tissue dictionary and a gene list (metadata)
def crawl(app):
	print("[i] reading unigene query from {0}".format(app["input_file"]))

	with open(app["input_file"]) as filehandler:
		# this list contains the metadata
		# for all queried genes from input_file
		genes = []
		
		# contains all fetched tissues		
		tissues = collections.Counter({})

		for line in filehandler:

			# sanitize input
			line = line.strip()

			# skip header line
			if line.startswith(">"):
				continue

			# build query
			query = build_query(line)

			# send request
			content, err = send_request(app, query)
			if err:
				print("[x] skipped, because of: {0}".format(content))
				continue

			# when query was with unigene id (ugid), parse cluster id (cid)
			# and organism (org) and send another request
			if "ugid" in query.keys():
				# check if there were no items found
				ident = re.compile(REGEX_ENTRYIDENT, re.S).findall(content)
		
				if not ident:
					print("[x] skipped, could not find an identifier line (REGEX_ENTRYIDENT)")
					continue

				identifier = ident[0].split("&nbsp;&nbsp;&nbsp;")
				if not identifier:
					print("[x] skipped, could not find identifier (&nbsp;&nbsp;&nbsp;)")
					continue

				cluster_information = identifier[1].strip().split(" ")
				organism, cluster_id = cluster_information[1].strip().split(".")

				# send another request
				query = build_query("[org]{0} [cid]{1}".format(organism, cluster_id))
				content, err = send_request(app, query)
				if err:
					print("[x] skipped, because of: {0}".format(content))
					continue

			# find tissues in response body
			crawled_tissues = parse_tissues(content)

			# add to existing tissues
			tissues = tissues + crawled_tissues

			# find meta data
			gene = parse_metadata(content)
			gene["tissue_cnt"] = reduce(lambda acc, curr: acc + curr, crawled_tissues.values(), 0)
			gene["query"] = query

			genes.append(gene)

		return genes, tissues

# Parses all tissues out of the given HTML
# and groups them by name
def parse_tissues(content):
	# find table with tissues
	table = re.compile(REGEX_TISSUETABLE, re.S).findall(content)

	if not table:
		print("[x] could not find a tissue table (REGEX_TISSUETABLE)")
		return collections.Counter({})

	table = table[0]
	tissues = re.compile(REGEX_TISSUE, re.X).findall(table)

	return collections.Counter(tissues)

# Parses genes / metadata out of the given HTML
# and returns a metadata object
# Meta data consists of ugid, cluster_id, organism, name, number of tissues, query line
def parse_metadata(content):
	# find identifier line
	ident = re.compile(REGEX_ENTRYIDENT, re.S).findall(content)
	if not ident:
		print("[x] could not find an identfier line (REGEX_ENTRYIDENT)")
		return {}

	identifier = ident[0].split("&nbsp;&nbsp;&nbsp;")
	if not identifier:
		print("[x] could not find valid identifier (&nbsp;&nbsp;&nbsp;)")
		return {}

	ugid = identifier[0].strip().split(":")[1]
	organism = identifier[1].strip().split(" ")[1].split(".")[0]
	cluster_id = identifier[1].strip().split(" ")[1].split(".")[1]

	# check if there were no items found
	name = re.compile(REGEX_ENTRYTITLE, re.S).findall(content)
	if not name:
		print("[x] could not find a title (REGEX_ENTRYTITLE)")
		return {}

	name = name[0]

	return {
		"ugid": ugid,
		"organism": organism,
		"cluster_id": cluster_id,
		"name": name,
	}

# Sends with given url query parameters a http request to
# the unigene database
def send_request(app, query):
	# prepare the request depending on query parameter
	if "link" in query.keys():
		
		# take the url query of the link and parse it to a dictionary
		# which will be passed to build_query
		query_dict = urllib.parse.parse_qs(urllib.parse.urlparse(query["link"]).query)

		# create query line
		query_line = ""
		for key, value in query_dict.items():
			# the query values are saves as array
			# as there might be more values to one key in the url query
			query_line = query_line + "[{0}]{1} ".format(key, value[0])

		query = build_query(query_line)

	if "org" in query.keys() and "cid" in query.keys():
		query["org"] = query["org"].lower().capitalize()
		query["maxest"] = 1000000000 if "maxest" not in query.keys() else int(query["maxest"])

	# encode query
	raw_q = urllib.parse.urlencode(query)
	q = raw_q.encode("utf-8")

	print("[i] crawling unigene with {0}".format(query))
	if app["debug"]:
		print("[d] URL: {0}?{1}".format(UNIGENE_URL, raw_q))
	try:
		# send HTTP GET request
		request = urllib.request.Request(url=UNIGENE_URL, data=q, method="POST")
		response = urllib.request.urlopen(request, timeout=app["timeout"])
		if response.status != 200:
			return "failed with {0}: {1}".format(response.status, response.reason), True

		content = response.read().decode("utf-8")

	except (urllib.error.HTTPError, urllib.error.URLError):
		return "could not connect to {0}".format(UNIGENE_URL), True

	except socket.timeout:
		return "could not connect to {0}: host does not respond".format(UNIGENE_URL), True

	# check if there were no items found
	title = re.compile(REGEX_ENTRYTITLE, re.S).findall(content)
	if not title:
		return "could not find an entry title (REGEX_ENTRYTITLE)", True

	if title[0].lower() == "no items found.":
		return "no unigene entry has been found", True

	if title[0].lower().endswith("has been retired"):
		return "this gene has been retired", True

	# return valid content
	return content, False

# Parses and build a query dictionary out of the given string line
def build_query(line):
	line = line.strip()
	args = line.split(" ")

	query = {}
	for arg in args:
		# split and filter out empty entries
		splitted = re.split(r"\[(.+?)\]", arg)
		splitted = list(filter(None, splitted))

		# skip empty lines
		if not splitted:
			continue

		# default lowercase
		key = splitted[0].lower().strip()
		value = splitted[1].lower().strip()

		# append to query dictionary
		query[key] = value

	# todo: validate
	return query

# Parses optional command line arguments
# and set the fields for the application
def parse_cmd_args(app):
	# program parameter passed by cmd
	args = sys.argv

	app["input_file"] = args[1] if len(args) > 1 else app["input_file"]
	app["output_dir"] = args[2] if len(args) > 2 else app["output_dir"]
	app["timeout"] = int(args[3]) if len(args) > 3 else app["timeout"]
	app["debug"] = True if len(args) > 4 and args[4] == "-debug" else app["debug"]

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

# Print all given tissues of all proteins to the console
def print_to_console(tissues, app):
	print()
	print("Tissues in {0}".format(app["name"]))
	print("="*(13 + len(app["name"])))

	for tissue, number in sorted(tissues.items()):
		print("{0}: {1}".format(tissue, number))

	print()

# Writes the found tissues of all ESTs to a .csv file.
def write_result_to_csv(tissues, app):
	path = os.path.join(app["path"], "tissues.csv")
	with open(path, "w") as filehandler:
		filehandler.write("Tissues in {0}:\n".format(app["name"]))
		filehandler.write("tissue;number\n")

		for tissue, number in sorted(tissues.items()):
			filehandler.write("{0};{1}\n".format(tissue, number))

	print("[i] csv saved to {0}".format(path))

# Writes the found metadata for each result to the .csv file.
# This meta data consists of ugid, cluster_id, organism, name, number of tissues, query line
def write_meta_to_csv(genes, app):
	path = os.path.join(app["path"], "meta.csv")
	with open(path, "w") as filehandler:
		filehandler.write("Meta data for {0}:\n".format(app["name"]))
		filehandler.write("ugid;cluster_id;organism;tissue_cnt;name;query\n")

		for gene in genes:
			if not gene:
				continue

			filehandler.write("{0};{1};{2};{3};{4};{5}\n".format(
				gene["ugid"], 
				gene["cluster_id"],
				gene["organism"],
				gene["tissue_cnt"],
				gene["name"],
				gene["query"]))

	print("[i] meta data saved to {0}".format(path))

def main():
	# define application, default program parameter
	app = {
		"input_file": "example",
		"output_dir": "output/",
		"timeout": 15,
		"debug": True,
	}

	# parse command line arguments
	parse_cmd_args(app)

	# validate command line arguments
	message, err = validate_cmd_args(app)
	if err:
		print("[x] {0}".format(message))
		sys.exit(1)

	# start crawling
	print("[i] starting unigene tissue crawler")
	genes, tissues = crawl(app)
	tissues = collections.OrderedDict(sorted(tissues.items()))

	# print the data to console
	print_to_console(tissues, app)

	# save meta data into a .csv file
	write_meta_to_csv(genes, app)

	# save result data (mutations) into a .csv file
	write_result_to_csv(tissues, app)

if __name__ == "__main__":
	main()