import os, sys

directory = sys.argv[1]
target_string = ""
if len(sys.argv) > 2:
	target_string = sys.argv[2]

if not os.path.isdir(directory):
	print("'%s' does not exist!"%directory)
	sys.exit(0)

files = sorted(os.listdir(directory))
if not files:
	print("'%s' is empty!" % directory)
	sys.exit(0)

for filename in files:
	full_filename = os.path.join(directory, filename)
	if os.stat(full_filename).st_size == 0:
		print("'%s' is empty!"%full_filename)
	else:
		with open(full_filename, "r") as f:
			print("'%s':"%full_filename)
			for line in f.readlines():
				if target_string in line:
					print("\t"+line.strip("\n"))
