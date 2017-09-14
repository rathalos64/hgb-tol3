# This scripts represents a deserializer for binary-saved matplotlib figures.
#
# Use this script to show the saved .matpltlib file as a figure 
# within your operating system.
# 
# In order to use this script, some python packages must be installed beforehand.
# Please refer to the "requirements.txt" one directory above.
# The best way to install the packages would be with pip:
# $ pip install -r requirements.txt
#
# Pass the name of the file simply as the first argument to the program
# Example:
# $ python showgraph.py ../output/nameofresult/example.matpltlib
#
# This script was developed with Python 3.6.2

# General imports
import sys
import os
import math

def main():
	args = sys.argv
	if len(args) == 1:
		print("[x] no arguments given. Abort...")
		sys.exit(1)

	path = args[1]
	if not os.path.exists(path):
		print("[x] input file {0} does not exist. Abort..".format(path))
		sys.exit(1)

	if os.stat(path).st_size == 0:
		print("[x] input file {0} is empty. Abort..".format(path))
		sys.exit(1)

	if not path.endswith(".matpltlib"):
		print("[x] input file {0} is not a valid .matpltlib file. Abort...".format(path))
		sys.exit(1)

	import matplotlib.pyplot as plt
	import dill

	with open(path, "rb") as filehandler:
		ax = dill.load(filehandler)

	plt.show()

if __name__ == "__main__":
	main()
