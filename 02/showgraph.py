#!/usr/bin/env python

import sys
import os
import math

def main():
	args = sys.argv
	if len(args) == 1:
		print("[x] no arguments given. Abort...")
		sys.exit(1)

	path = args[1]
	if os.stat(path).st_size == 0:
		print("[x] path {0} is either empty or does not exist. Abort..".format(path))
		sys.exit(1)

	if not path.endswith(".matpltlib"):
		print("[x] path {0} not a valid .matpltlib file. Abort...".format(path))
		sys.exit(1)

	import matplotlib.pyplot as plt
	import dill

	with open(path, "rb") as filehandler:
		ax = dill.load(filehandler)

	plt.show()

if __name__ == "__main__":
	main()
