import os
import sys

import re


def find_floats_in_string(string):
	return list(map(float, re.findall(r'\d +\.\d +', string)))

def get_avg_runtime(path, filename):
	output = open(os.path.join(path, filename), "r").read()
	if not output:
		return None
	runtimes = find_floats_in_string(output) 
	return sum(runtimes)/len(runtimes)

def outputs2csv():
	results = {}
	for filename in os.listdir(OUTPUTS_PATH):
		if filename.split(".")[-1] != "dat":
			print("incompatible file: %s" % filename)
			continue
		print("extracting data from %s" % filename)
		results[filename.split(".")[0]] = get_avg_runtime(OUTPUTS_PATH, filename)

	print(results)
	#with open(os.path.join(DATA_PATH, filename.replace(".dat", ".csv")), "w") as f:

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("wrong number of arguments: "+str(sys.argv[1:]))
		sys.exit(-1)

	OUTPUTS_PATH = sys.argv[1]
	DATA_PATH = sys.argv[2]

	outputs2csv()
