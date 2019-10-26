### PROGRAMS ###
SEQUENTIAL = seq_1D
PARALLEL = mpi_1D
LANGUAGE = c
#LANGUAGE = cpp

GRAPH = graph.py
PRINT_FILES = print_files.py


### PARAMETERS for multiple mpi runs ###
# --> use 'make run' for submitting multiple jobs with different configurations
NUM_RUNS = 1

#program will be run with various mpi configurations each specified by #ranks, slot distribution across nodes/cpus and problem sizes
NUM_RANKS_ARRAY = 4 8
#NUM_RANKS_ARRAY = 1 2 4 8 16
SLOT_DISTRIBUTE_ARRAY = fillup
#SLOT_DISTRIBUTE_ARRAY = fillup 1perhost 2perhost
PROBLEM_SIZE_ARRAY = 2000 4000

### PARAMETERS for single mpi run ###
# --> use 'make run_mpi' for debugging with these default values
#problem_size := $(word 1, $(PROBLEM_SIZE_ARRAY))
#num_ranks := $(word 1, $(NUM_RANKS_ARRAY))
#slot_distribute := $(word 1, $(SLOT_DISTRIBUTE_ARRAY))


### OUTPUT + ERRORS ###
OUTPUTS_DIR = ./outputs
ERRORS_DIR = ./errors
RESULTS_DIR = ./results

SEQ_NAME = $(SEQUENTIAL)_$${problem_size}
MPI_NAME = $(PARALLEL)_$${problem_size}_$${slot_distribute}_$${num_ranks}_slots

SEQ_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(SEQ_NAME).dat
MPI_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(MPI_NAME).dat

#create errors in cwd to easily spot them!
SEQ_ERROR_FILENAME = $(ERRORS_DIR)/$(SEQ_NAME).err
MPI_ERROR_FILENAME = $(ERRORS_DIR)/$(MPI_NAME).err

PYTHON = python
GRAPH = graph.py


### COMPILERS + FLAGS ###
CC = gcc
CC_FLAGS = -O2 -std=c99 -Wall

MPI = mpicc									#mpi++ for c++!
MPI_FLAGS = -O2 -std=c99 -Wall


### Windows/LCC2 SPECIFIC SETUP ###
ifeq ($(OS), Windows_NT)
	MPI_VERSION = 3.1.3

	# additional linking required on Windows
	#@Marc: add Windows environmental variable 'CYGWIN_PATH --> C:\cygwin64'
	MPI_FLAGS += -I "$(CYGWIN_PATH)\usr\include" 
	MPI_FLAGS += -I "$(CYGWIN_PATH)\lib\openmpi"

	SEQ_EXEC_COMMAND = echo "running '$(SEQ_NAME)'";
	SEQ_EXEC_COMMAND += ./$(SEQUENTIAL) $${problem_size} 1>> $(SEQ_OUTPUT_FILENAME) 2>> $(SEQ_ERROR_FILENAME)

	MPI_EXEC_COMMAND = echo "running '$(MPI_NAME)'";
	MPI_EXEC_COMMAND += mpiexec -n $${num_ranks} --oversubscribe ./$(PARALLEL) $${problem_size} 1>> $(MPI_OUTPUT_FILENAME) 2>> $(MPI_ERROR_FILENAME)
else
	MPI_VERSION = 4.0.1

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
	SGE_MPI_FLAGS += -pe openmpi-$(slot_distribute) $(num_ranks)
	#cannot use '-b yes' because openmpi module needs to be loaded on lcc2 cluster


	SEQ_EXEC_COMMAND = qsub $(SGE_SEQ_FLAGS) ./$(SEQUENTIAL) $${problem_size}

	MPI_EXEC_COMMAND = echo -e "#!/bin/bash\nmodule load openmpi/$(MPI_VERSION) 2>/dev/null\nmpiexec -n $${num_ranks} ./$(PARALLEL) $${problem_size}" > job.tmp
endif

### CREATE NON-EXISTING DIRECTORIES ###
$(OUTPUTS_DIR):
	@mkdir -p $@

$(ERRORS_DIR):
	@mkdir -p $@

$(RESULTS_DIR):
	@mkdir -p $@


### COMPILE ###
.PHONEY: compile_seq
compile_seq: $(SEQUENTIAL).$(LANGUAGE)
	@$(CC) $(CC_FLAGS) $< -o $(SEQUENTIAL)

.PHONEY: compile_mpi
compile_mpi: $(PARALLEL).$(LANGUAGE)
ifeq ($(OS), Windows_NT)
	@$(MPI) $(MPI_FLAGS) $< -o $(PARALLEL)
else
	@module load openmpi/$(MPI_VERSION) 2>/dev/null && $(MPI) $(MPI_FLAGS) $< -o $(PARALLEL)
endif

.PHONEY: all
all: compile_seq compile_mpi


### RUN ###
.PHONEY: run_seq run_seq_single
run_seq_single: PROBLEM_SIZE_ARRAY := $(word 1, $(PROBLEM_SIZE_ARRAY)) 

run_seq run_seq_single: compile_seq $(OUTPUTS_DIR) $(ERRORS_DIR)
	@for problem_size in $(PROBLEM_SIZE_ARRAY); do \
		$(SEQ_EXEC_COMMAND); \
	done
	@echo


.PHONEY: run_mpi run_mpi_single
run_mpi_single: PROBLEM_SIZE_ARRAY := $(word 1, $(PROBLEM_SIZE_ARRAY)) 
run_mpi_single: NUM_RANKS_ARRAY := $(word 1, $(NUM_RANKS_ARRAY)) 
run_mpi_single: SLOT_DIRSTRIBUTE_ARRAY := $(word 1, $(SLOT_DISTRIBUTE_ARRAY))

run_mpi run_mpi_single: compile_mpi $(OUTPUTS_DIR) $(ERRORS_DIR)
	for problem_size in $(PROBLEM_SIZE_ARRAY); do \
		for num_ranks in $(NUM_RANKS_ARRAY); do \
			for slot_distribute in $(SLOT_DISTRIBUTE_ARRAY); do \
				$(MPI_EXEC_COMMAND); \
			done \
		done \
	done
	@echo


.PHONEY: run
run: super_clean run_seq run_mpi check_errors check_outputs

### VISUALIZE ###
# TODO --> runtime, speedup, efficiency, ...
.PHONEY: graphs
graphs: $(GRAPH) $(OUTPUTS_DIR) $(RESULTS_DIR)
	@echo creating some fancy plots...
	@$(PYTHON) ./$(GRAPH) $(OUTPUTS_DIR) $(RESULTS_DIR)
	@echo 

### CHECK OUTPUTS & ERRORS ###
.PHONEY: check_outputs
check_outputs:
	@$(PYTHON) ./$(PRINT_FILES) $(OUTPUTS_DIR) seconds
	@echo 

.PHONEY: check_errors
check_errors:
	@$(PYTHON) ./$(PRINT_FILES) $(ERRORS_DIR)
	@echo 

### CLEANUP ###
# '|| true' --> ignore errors if files don't exist
.PHONEY: clean
clean:
	@rm -r $(ERRORS_DIR) 2>/dev/null  || true
ifeq ($(OS), Windows_NT)
	@rm $(SEQUENTIAL).exe $(PARALLEL).exe 2>/dev/null || true
else
	@rm $(SEQUENTIAL) $(PARALLEL) 2>/dev/null  || true
endif
	@echo "your directory has been cleaned of executables and errors!"


.PHONEY: super_clean
super_clean: clean
	@rm -r $(OUTPUTS_DIR) 2>/dev/null  || true
	@echo "... and also outputs!"
	@echo 
