import os
import sys

import re
import csv

def dict2csv(dim_dict, filename):
	with open(os.path.join(DATA_PATH, filename), 'w', newline='') as csvfile:
		fieldnames = ["dimensions"]
		problem_sizes = []
		for col_name in dim_dict:
			fieldnames.append(col_name)
			for ps in dim_dict[col_name]:
				if ps not in problem_sizes:
					problem_sizes.append(ps)

		problem_sizes = sorted(problem_sizes)
				
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()

		for ps in problem_sizes:
			row = {"dimensions": ps}
			for col_name in dim_dict:
				#print(col_name, dim_dict[col_name])
				if dim_dict[col_name] is not None and ps in dim_dict[col_name]:
					row[col_name] = dim_dict[col_name][ps]
				else:
					dim_dict[col_name] = None
			
			writer.writerow(row)


def find_ints_in_string(string):
	return list(map(int, re.findall(r'\d+', string)))

def find_floats_in_string(string):
	return list(map(float, re.findall(r'\d+\.\d+', string)))

def get_avg_runtime(path, filename):
	output = open(os.path.join(path, filename), "r").read()
	if not output:
		return None

	str_runtimes = " - ".join([s.split("seconds")[0] for s in output.split("The process took") if "seconds" in s])
	runtimes = find_floats_in_string(str_runtimes)
	if len(runtimes) == 0:
		print("missing value for '%s'"%filename)
		return None

	return sum(runtimes)/len(runtimes)

def outputs2csv():
	results = {}
	for filename in os.listdir(OUTPUTS_PATH):
		if filename.split(".")[-1] != "dat":
			print("incompatible file: %s" % filename)
			continue

		# categorise different output filenames into groups each comprising of average runtimes for different problem_sizes!
		basename = filename.split(".")[0]
		if basename.count("_") == 2:
			#seq category
			group_name, problem_size = basename.rsplit("_", 1)
		elif basename.count("_") in [3, 4]:
			group_name, problem_size, slot_distribute, num_slots = basename.rsplit("_", 3)
			group_name = "_".join([group_name, num_slots])
		else:
			print("unsupported_basename: %s!\n"%basename)
			sys.exit(1)

		problem_size = int(problem_size)
		if group_name not in results:
			results[group_name] = {}
		#print("extracting data from %s" % filename)

		results[group_name][problem_size] = get_avg_runtime(OUTPUTS_PATH, filename)

	#print(results)
	dict2csv(results, "Real.csv")
	

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("wrong number of arguments: "+str(sys.argv[1:]))
		sys.exit(-1)

	OUTPUTS_PATH = sys.argv[1]
	DATA_PATH = sys.argv[2]

	outputs2csv()
