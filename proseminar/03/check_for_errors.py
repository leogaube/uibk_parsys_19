import os, sys

error_dir = sys.argv[1]

if not os.path.isdir(error_dir):
	print("'%s' does not exist!"%error_dir)
	sys.exit(0)

error_files = os.listdir(error_dir)
if not error_files:
	print("'%s' is empty :)" % error_dir)
	sys.exit(0)

for filename in error_files:
	if filename[-4:] != ".err":
		continue

	full_filename = os.path.join(error_dir, filename)
	if os.stat(full_filename).st_size == 0:
		print("'%s' is empty :)"%full_filename)
	else:
		with open(full_filename, "r") as f:
			print("'%s' contains errors:"%full_filename)
			for line in f.readlines():
				print("\t"+line)
