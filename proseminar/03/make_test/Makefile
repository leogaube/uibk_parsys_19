### PROGRAMS ###
SEQUENTIAL_ARRAY = seq_3D seq_2D seq_1D
PARALLEL_ARRAY = mpi_1D mpi_2D mpi_3D
DEPENDENCY_ARRAY = heat_stencil verification
DEFINES = -DVERBOSE -DSLICE_3D
LANGUAGE = c
#LANGUAGE = cpp

PREPROCESS = preprocess.py
PLOT = plot.py
PRINT_FILES = print_files.py

### CONFIGURATIONS ###
# program will be run with various mpi configurations 
# each specified by #ranks, slot distribution across nodes/cpus and problem sizes 
# --> use 'make run' for submitting multiple jobs with different configurations
# the default configuration consists of all first elements of the *_ARRAY variables
# --> use 'make run_mpi_single' for submitting a single job. 
NUM_RUNS = 1

NUM_RANKS_ARRAY = 4 8
#NUM_RANKS_ARRAY = 1 2 4 8 16
SLOT_DISTRIBUTE_ARRAY = fillup
#SLOT_DISTRIBUTE_ARRAY = fillup 1perhost 2perhost
PROBLEM_SIZE_ARRAY = 8 16 32


### OUTPUT + ERRORS ###
OUTPUTS_DIR = ./outputs
ERRORS_DIR = ./errors
DATA_DIR = ./data
RESULTS_DIR = ./results

SEQ_NAME = $${sequential}_$${problem_size}
MPI_NAME = $${parallel}_$${problem_size}_$${slot_distribute}_$${num_ranks}_slots

SEQ_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(SEQ_NAME).dat
MPI_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(MPI_NAME).dat

SEQ_ERROR_FILENAME = $(ERRORS_DIR)/$(SEQ_NAME).err
MPI_ERROR_FILENAME = $(ERRORS_DIR)/$(MPI_NAME).err


######## SETUP ##########
### COMPILERS + FLAGS ###
CC = gcc
CC_FLAGS = -O2 -std=c99 -Wall

MPI = mpicc									#mpi++ for c++!
MPI_FLAGS = -O2 -std=c99 -Wall


PYTHON = python

### Windows/LCC2 specific setup ###
ifeq ($(OS), Windows_NT)
	MPI_VERSION = 3.1.3

	LOAD_MPI = ""

	DEP_OBJECTS = $(addsuffix .o, $(DEPENDENCY_ARRAY))
	SEQ_EXECUTABLES = $(addsuffix .exe, $(SEQUENTIAL_ARRAY))
	MPI_EXECUTABLES = $(addsuffix .exe, $(PARALLEL_ARRAY))
	EXECUTABLE_PATTERN = %.exe

	# additional linking required on Windows
	# Windows environmental variable 'CYGWIN_PATH --> C:\cygwin' required
	MPI_FLAGS += -I "$(CYGWIN_PATH)\usr\include" 
	MPI_FLAGS += -I "$(CYGWIN_PATH)\lib\openmpi"

	#[ $${sequential} = seq_1D ] && problem_size="$$(($${problem_size}*$${problem_size}*$${problem_size}))"; \
	[ $${sequential} = seq_2D ] && problem_size="$$(($${problem_size}*$${problem_size}))"; \

	SEQ_EXEC_COMMAND = echo "running '$(SEQ_NAME)'";
	SEQ_EXEC_COMMAND += ./$${sequential} $${problem_size} 1>> $(SEQ_OUTPUT_FILENAME) 2>> $(SEQ_ERROR_FILENAME)

	MPI_EXEC_COMMAND = [ $${parallel} = mpi_1D ] && problem_size=eval "awk -v num=$${problem_size} 'BEGIN{print int(num^(3))}'";
	MPI_EXEC_COMMAND += [ $${parallel} = mpi_2D ] && problem_size=eval "awk -v num=$${problem_size} 'BEGIN{print int(num^(3/2))}'";
	MPI_EXEC_COMMAND += echo "running '$(MPI_NAME)'";
	MPI_EXEC_COMMAND += mpiexec -n $${num_ranks} --oversubscribe ./$${parallel} $${problem_size} 1>> $(MPI_OUTPUT_FILENAME) 2>> $(MPI_ERROR_FILENAME)
else
	MPI_VERSION = 4.0.1

	LOAD_MPI = module load openmpi/$(MPI_VERSION) 2>/dev/null && 

	DEP_OBJECTS = $(addsuffix .o, $(DEPENDENCY_ARRAY))
	SEQ_EXECUTABLES = $(SEQUENTIAL_ARRAY)
	MPI_EXECUTABLES = $(PARALLEL_ARRAY)
	EXECUTABLE_PATTERN = %

	#TODO add and test with '-l h_rt= / h_vmem=' parameters
	SGE_SEQ_FLAGS = -q std.q
	SGE_SEQ_FLAGS += -cwd
	SGE_SEQ_FLAGS += -N $(SEQ_NAME)
	SGE_SEQ_FLAGS += -o $(SEQ_OUTPUT_FILENAME)
	SGE_SEQ_FLAGS += -e $(SEQ_ERROR_FILENAME)
	SGE_SEQ_FLAGS += -b yes

	SGE_MPI_FLAGS = -q std.q
	SGE_MPI_FLAGS += -cwd
	SGE_MPI_FLAGS += -N $(MPI_NAME)
	SGE_MPI_FLAGS += -o $(MPI_OUTPUT_FILENAME)
	SGE_MPI_FLAGS += -e $(MPI_ERROR_FILENAME)
	SGE_MPI_FLAGS += -pe openmpi-$$slot_distribute $$num_ranks
	#cannot use '-b yes' because openmpi module needs to be loaded on lcc2 cluster TODO

	SEQ_EXEC_COMMAND = [ $${sequential} = seq_1D ] && problem_size=eval "awk -v num=$${problem_size} 'BEGIN{print int(num^(3))}'";
	SEQ_EXEC_COMMAND += [ $${sequential} = seq_2D ] && problem_size=eval "awk -v num=$${problem_size} 'BEGIN{print int(num^(3/2))}'";
	SEQ_EXEC_COMMAND += echo $${problem_size};
	SEQ_EXEC_COMMAND += qsub $(SGE_SEQ_FLAGS) ./$${sequential} $${problem_size}

	MPI_EXEC_COMMAND = [ $${parallel} = mpi_1D ] && problem_size=eval "awk -v num=$${problem_size} 'BEGIN{print int(num^(3))}'";
	MPI_EXEC_COMMAND += [ $${parallel} = mpi_2D ] && problem_size=eval "awk -v num=$${problem_size} 'BEGIN{print int(num^(3/2))}'";
	MPI_EXEC_COMMAND += echo -e "module load openmpi/$(MPI_VERSION) 2>/dev/null\nmpiexec -n $${num_ranks} ./$${parallel} $${problem_size}" > job.tmp;
	MPI_EXEC_COMMAND += qsub $(SGE_MPI_FLAGS) job.tmp;
	MPI_EXEC_COMMAND += rm job.tmp
endif


############## TARGETS ################
.DEFAULT_GOAL := run_mpi


.PHONEY: run_mpi run_mpi_single
run_mpi_single: PARALLEL_ARRAY := $(word 1, $(PARALLEL_ARRAY)) 
run_mpi_single: PROBLEM_SIZE_ARRAY := $(word 1, $(PROBLEM_SIZE_ARRAY)) 
run_mpi_single: NUM_RANKS_ARRAY := $(word 1, $(NUM_RANKS_ARRAY)) 
run_mpi_single: SLOT_DIRSTRIBUTE_ARRAY := $(word 1, $(SLOT_DISTRIBUTE_ARRAY))

run_mpi run_mpi_single:
	for parallel in $(PARALLEL_ARRAY); do \
		for run in in $(seq 1 $(NUM_RUNS)); do \
			for problem_size in $(PROBLEM_SIZE_ARRAY); do \
				for num_ranks in $(NUM_RANKS_ARRAY); do \
					for slot_distribute in $(SLOT_DISTRIBUTE_ARRAY); do \
						$(MPI_EXEC_COMMAND); \
					done \
				done \
			done \
		done \
	done
	@echo