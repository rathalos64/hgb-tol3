#!/usr/bin/env python

from Bio import Entrez
from Bio import SeqIO

import json

Entrez.email = "s1510458019@students.fh-hagenberg.at"

class Terminal:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def main():
	objects = []
	
	# crawl all genes from given input
	with open('input') as fh:
		for line in fh.readlines():
			objects.append(crawl(line))

	# plot single statistic for each found gene
	for obj in objects:
		plot_single(obj)

def crawl(id):
	# query object
	obj = {}
	obj['id'] = id
	obj['polyA'] = {}
	obj['polyA']['location'] = {}
	obj['polyA']['distribution'] = {}

	# query entrez
	with Entrez.efetch(db='nucleotide', id=obj['id'], rettype='gb', retmode='xml') as handler:
		record = Entrez.read(handler)[0]
		obj['name'] = record['GBSeq_definition']
		
		# retrieve base sequence
		obj['sequence'] = record['GBSeq_sequence']

		# look up for regulatory feature and get location
		feature = list(filter(lambda x: x['GBFeature_key'] == 'regulatory', record['GBSeq_feature-table']))[0]
		location = feature['GBFeature_location'].split('..')
		obj['polyA']['location']['start'] = int(location[0])
		obj['polyA']['location']['end'] = int(location[1])

		# look up for polyA-site feature and get distance to tail
		feature = list(filter(lambda x: x['GBFeature_key'] == 'polyA_site', record['GBSeq_feature-table']))[0]
		obj['polyA']['distance'] = int(feature['GBFeature_location']) - obj['polyA']['location']['end']

		# get taxon number
		feature = list(filter(lambda x: x['GBFeature_key'] == 'source', record['GBSeq_feature-table']))[0]
		qualifier = list(filter(lambda x: x['GBQualifier_name'] == 'db_xref', feature['GBFeature_quals']))[0]
		obj['family'] = qualifier['GBQualifier_value']

		handler.close()

	# base polyA data
	obj['polyA']['sequence'] = obj['sequence'][obj['polyA']['location']['start'] - 1:obj['polyA']['location']['end'] - 1]
	obj['polyA']['length'] = len(obj['polyA']['sequence'])

	# percentual polyA position
	obj['polyA']['position'] = round(obj['polyA']['location']['start'] / len(obj['sequence']), 2)

	# printable polyA
	obj['polyA']['print'] = obj['sequence'][obj['polyA']['location']['start'] - 11:obj['polyA']['location']['start'] - 1]
	obj['polyA']['print'] += Terminal.BOLD + obj['polyA']['sequence'].upper() + Terminal.END
	obj['polyA']['print'] += obj['sequence'][obj['polyA']['location']['end'] - 1:obj['polyA']['location']['end'] + 10]

	# calculate distribution
	obj['polyA']['distribution']['C'] = obj['polyA']['sequence'].count('c') / obj['polyA']['length']
	obj['polyA']['distribution']['G'] = obj['polyA']['sequence'].count('g') / obj['polyA']['length']
	obj['polyA']['distribution']['T'] = obj['polyA']['sequence'].count('t') / obj['polyA']['length']
	obj['polyA']['distribution']['A'] = obj['polyA']['sequence'].count('a') / obj['polyA']['length']

	return obj

def plot_single(obj):
	print(f"####### [{obj['id']}]: {obj['name']} #######################################")
	print()

	print(f"family: {obj['family']}")
	print()
	print(f"retrieved polyA-Signal-Sequence: {obj['polyA']['sequence'].upper()}")
	print(f"len: {obj['polyA']['length']}")
	print(f"location: [{obj['polyA']['location']['start']}..{obj['polyA']['location']['end']}]")
	print(f"in sequence: ...{obj['polyA']['print']}...")
	print()

	print(f"percentual position in sequence: {obj['polyA']['position'] * 100}%")
	print()

	print("distribution:")
	print(f"\tC: {obj['polyA']['distribution']['C'] * 100}%")
	print(f"\tG: {obj['polyA']['distribution']['G'] * 100}%")
	print(f"\tT: {obj['polyA']['distribution']['T'] * 100}%")
	print(f"\tA: {obj['polyA']['distribution']['A'] * 100}%")
	
if __name__ == "__main__":
	main()