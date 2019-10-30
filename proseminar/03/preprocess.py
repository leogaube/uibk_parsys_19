import os
import sys

import re
import csv

def dict2csv(dim_dict, filename):
	with open(os.path.join(DATA_PATH, filename), 'w', newline='') as csvfile:
		fieldnames = ["room_size"]
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
			row = {"room_size": ps}
			for col_name in dim_dict:
				if ps in dim_dict[col_name]:
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
	if "seq" in filename:
		print(filename, str_runtimes)
	runtimes = find_floats_in_string(str_runtimes)
	if len(runtimes) == 0:
		return None

	return sum(runtimes)/len(runtimes)

def outputs2csv():
	results = {}
	for filename in os.listdir(OUTPUTS_PATH):
		if filename.split(".")[-1] != "dat":
			print("incompatible file: %s" % filename)
			continue

		basename = filename.split(".")[0]
		if basename.count("_") == 2:
			program_type, dim, problem_size = basename.split("_")
			name = program_type
		elif basename.count("_") == 4:
			program_type, dim, problem_size, slot_distribute, num_slots = basename.split("_")
			name = "_".join([program_type, num_slots, slot_distribute])
		else:
			print("unsupported_basename: %s"%basename)


		problem_size = int(problem_size)
		if dim not in results:
			results[dim] = {}
		if name not in results[dim]:
			results[dim][name] = {}
		print("extracting data from %s" % filename)

		results[dim][name][problem_size] = get_avg_runtime(OUTPUTS_PATH, filename)

	print(results)
	for dim in results:
		dict2csv(results[dim], dim+".csv")
	

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("wrong number of arguments: "+str(sys.argv[1:]))
		sys.exit(-1)

	OUTPUTS_PATH = sys.argv[1]
	DATA_PATH = sys.argv[2]

	outputs2csv()
