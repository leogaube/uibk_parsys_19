
SEQUENTIAL = heat_stencil_1D_seq
PARALLEL = heat_stencil_1D_mpi
LANGUAGE = c
#LANGUAGE = cpp


#multiple mpi runs with different configurations
NUM_RUNS = 1
PROBLEM_SIZE_START = 2000
PROBLEM_SIZE_END = 16000

NUM_RANKS_ARRAY = 4 2
#NUM_RANKS_ARRAY = 1 2 4 8 16
RANK_CONFIG_ARRAY = fillup
#RANK_CONFIG_ARRAY = fillup 1perhost 2perhost

#default values for single mpi run
PROBLEM_SIZE = $(PROBLEM_SIZE_START)
NUM_RANKS := $(word 1, $(NUM_RANKS_ARRAY))
RANK_CONFIG := $(word 1, $(RANK_CONFIG_ARRAY))

### SETUP ###

CYGWIN_PATH = C:\cygwin						#@Marc --> C:\cygwin64

CC = gcc
CC_FLAGS = -O2 -std=c99 -Wall

MPI = mpicc									#mpi++ for c++!
MPI_FLAGS = -O2 -std=c99 -Wall

#PYTHON=python


OUTPUTS_DIR = ./outputs
#RESULTS_DIR = ./results
SEQ_NAME = seq
MPI_NAME = $(NUM_RANKS)_$(RANK_CONFIG)

SEQ_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(SEQ_NAME).dat
MPI_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(MPI_NAME).dat

SEQ_ERROR_FILENAME = $(OUTPUTS_DIR)/$(SEQ_NAME).err
MPI_ERROR_FILENAME = $(OUTPUTS_DIR)/$(MPI_NAME).err

# Windows/LCC2 specific setup
ifeq ($(OS), Windows_NT)
	MPI_VERSION=3.1.3

	# additional linking required on Windows
	MPI_FLAGS += -I "$(CYGWIN_PATH)\usr\include" 
	MPI_FLAGS += -I "$(CYGWIN_PATH)\lib\openmpi"

	EXTENSION = .exe
else
	MPI_VERSION = 4.0.1

	SGE_SEQ_FLAGS = -q std.q
	SGE_SEQ_FLAGS += -cwd
	SGE_SEQ_FLAGS += -N $(SEQ_NAME)
	SGE_SEQ_FLAGS += -o $(SEQ_OUTPUT_FILENAME)
	SGE_SEQ_FLAGS += -e $(SEQ_ERROR_FILENAME)
	SGE_SEQ_FLAGS += -b

	SGE_MPI_FLAGS = -q std.q
	SGE_MPI_FLAGS += -cwd
	SGE_MPI_FLAGS += -N $(MPI_NAME)
	SGE_MPI_FLAGS += -o $(MPI_OUTPUT_FILENAME)
	SGE_MPI_FLAGS += -e $(MPI_ERROR_FILENAME)
	SGE_MPI_FLAGS += -pe openmpi-$(MPI_SETUP_TYPE)
	SGE_MPI_FLAGS += -b

	EXTENSION=.o
endif

# create non-existing directories
$(OUTPUTS_DIR):
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
	qsub $(SGE_MPI_FLAGS) mpiexec -n $(NUM_RANKS) ./$(PARALLEL) $(PROBLEM_SIZE)
endif

#TODO run mutiple times with multiple rank setups: (for \ ...)
.PHONEY: run_all_mpi
run_all_mpi: run_mpi

.PHONEY: run_all
run_all: run_seq run_all_mpi


### VISUALIZE ###

#.PHONEY graphs
#graphs: graph.py
#	@python  $(OUTPUTS_DIR) $(OUTPUTS_DIR)

### CLEANUP ###
.PHONEY: clean
clean:
	@rm *$(EXTENSION) || true			#hide output with: 2> /dev/null
	@rm -r $(OUTPUTS_DIR) || true
