### C/C++ SEQ/MPI PROGRAMS ###
SEQUENTIAL = seq_1D
PARALLEL = mpi_1D
#LANGUAGE = cpp
LANGUAGE = c


### PARAMETERS for multiple mpi runs ###
# --> 'make run_all' for submitting multiple jobs with different configurations
NUM_RUNS = 1

#problem size will be doubled for each consecutive run --> 2000 4000 8000
PROBLEM_SIZE_START = 2000
NEXT_PROBLEM_SIZE = 2*$(PROBLEM_SIZE_START)
NUM_PROBLEM_SIZES = 3

#program will be run with various mpi configurations each specified by #ranks and slot distribution across nodes/cpus
NUM_RANKS_ARRAY = 4
#NUM_RANKS_ARRAY = 1 2 4 8 16
SLOT_DISTRIBUTE_ARRAY = fillup
#SLOT_DISTRIBUTE_ARRAY = fillup 1perhost 2perhost


### PARAMETERS for single mpi run ###
# --> 'make run_mpi' for debugging with these default values
PROBLEM_SIZE := $(PROBLEM_SIZE_START)
NUM_RANKS := $(word 1, $(NUM_RANKS_ARRAY))
SLOT_DISTRIBUTE := $(word 1, $(SLOT_DISTRIBUTE_ARRAY))


### OUTPUT + ERRORS ###
OUTPUTS_DIR = ./outputs
RESULTS_DIR = ./results

SEQ_NAME = $(SEQUENTIAL)
MPI_NAME = $(PARALLEL)_$(SLOT_DISTRIBUTE)_$(NUM_RANKS)

SEQ_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(SEQ_NAME).dat
MPI_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(MPI_NAME).dat

SEQ_ERROR_FILENAME = $(OUTPUTS_DIR)/$(SEQ_NAME).err
MPI_ERROR_FILENAME = $(OUTPUTS_DIR)/$(MPI_NAME).err

PYTHON = python
GRAPH = graph.py


### COMPILERS + FLAGS ###
CC = gcc
CC_FLAGS = -O2 -std=c99 -Wall

MPI = mpicc									#mpi++ for c++!
MPI_FLAGS = -O2 -std=c99 -Wall


### Windows/LCC2 SPECIFIC SETUP ###
ifeq ($(OS), Windows_NT)
	MPI_VERSION=3.1.3

	# additional linking required on Windows
	#@Marc: add Windows environmental variable 'CYGWIN_PATH --> C:\cygwin64'
	MPI_FLAGS += -I "$(CYGWIN_PATH)\usr\include" 
	MPI_FLAGS += -I "$(CYGWIN_PATH)\lib\openmpi"
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
	SGE_MPI_FLAGS += -pe openmpi-$(SLOT_DISTRIBUTE) $(NUM_RANKS)
	SGE_MPI_FLAGS += -b yes
endif

### CREATE NON-EXISTING DIRECTORIES ###
.PHONEY: $(OUTPUTS_DIR)
$(OUTPUTS_DIR):
	@mkdir -p $@

.PHONEY: $(RESULTS_DIR)
$(RESULTS_DIR):
	@mkdir -p $@


### COMPILE ###
.PHONEY: compile_seq
compile_seq: $(SEQUENTIAL).$(LANGUAGE)
	@$(CC) $(CC_FLAGS) $< -o $(SEQUENTIAL)

.PHONEY: compile_mpi
compile_mpi: $(PARALLEL).$(LANGUAGE)
	@$(MPI) $(MPI_FLAGS) $< -o $(PARALLEL)

.PHONEY: all
all: compile_seq compile_mpi


### RUN ###
.PHONEY: run_seq
run_seq: compile_seq $(OUTPUTS_DIR)
ifeq ($(OS), Windows_NT)
	./$(SEQUENTIAL) $(PROBLEM_SIZE) 1>> $(SEQ_OUTPUT_FILENAME) 2>> $(SEQ_ERROR_FILENAME)
else
	qsub $(SGE_SEQ_FLAGS) ./$(SEQUENTIAL) $(PROBLEM_SIZE)
endif

.PHONEY: run_mpi
run_mpi: compile_mpi $(OUTPUTS_DIR)
ifeq ($(OS), Windows_NT)
	mpiexec -n $(NUM_RANKS) ./$(PARALLEL) $(PROBLEM_SIZE) 1>> $(MPI_OUTPUT_FILENAME) 2>> $(MPI_ERROR_FILENAME)
else
	qsub $(SGE_MPI_FLAGS) module load openmpi/4.0.1 && mpiexec -n $(NUM_RANKS) ./$(PARALLEL) $(PROBLEM_SIZE)
endif

#TODO run mutiple times with multiple rank setups: (for \ ...)
.PHONEY: run_all_mpi
run_all_mpi: run_mpi

.PHONEY: run_all
run_all: run_seq run_all_mpi


### VISUALIZE ###
# TODO --> runtime, speedup, efficiency, ...
.PHONEY: graphs
graphs: $(GRAPH) $(RESULTS_DIR)
	@$(PYTHON) ./$(GRAPH) $(OUTPUTS_DIR) $(RESULTS_DIR)


### CLEANUP ###
# ignore errors if files don't exist: use '|| true' 
.PHONEY: clean
clean:
	@rm -r $(OUTPUTS_DIR) || true
ifeq ($(OS), Windows_NT)
	@rm *.exe || true
else
	# all executable files!
	@ls | grep -v "\." | xargs rm || true
endif
