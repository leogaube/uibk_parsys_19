### PROGRAMS ###
SEQUENTIAL_ARRAY = seq_3D
PARALLEL_ARRAY = cubes
DEPENDENCY_ARRAY = heat_stencil
DEFINES = -D VERBOSE -D SLICE_3D
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

NUM_RANKS_ARRAY = 4 8 16 2 1
#NUM_RANKS_ARRAY = 1 2 4 8 16
SLOT_DISTRIBUTE_ARRAY = fillup
#SLOT_DISTRIBUTE_ARRAY = fillup 1perhost 2perhost
PROBLEM_SIZE_ARRAY = 32 64 128


### OUTPUT + ERRORS ###
OUTPUTS_DIR = ./outputs
ERRORS_DIR = ./errors
DATA_DIR = ./data
RESULTS_DIR = ./results

SEQ_NAME = $${sequential}_$${problem_size}
MPI_NAME = $${parallel}_$${problem_size}_$${slot_distribute}_$${num_ranks}-slots

SEQ_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(SEQ_NAME).dat
MPI_OUTPUT_FILENAME = $(OUTPUTS_DIR)/$(MPI_NAME).dat

SEQ_ERROR_FILENAME = $(ERRORS_DIR)/$(SEQ_NAME).err
MPI_ERROR_FILENAME = $(ERRORS_DIR)/$(MPI_NAME).err


######## SETUP ##########
### COMPILERS + FLAGS ###
CC = gcc
CC_FLAGS = -O2 -std=c99 -Wall -g

MPI = mpicc									#mpi++ for c++!
MPI_FLAGS = -O2 -std=c99 -Wall -g
LIBRARIES = -lm

PYTHON = python

### Executables ###
DEP_OBJECTS = $(addsuffix .o, $(DEPENDENCY_ARRAY))
SEQ_EXECUTABLES = $(addsuffix $(EXECUTABLE_EXT), $(SEQUENTIAL_ARRAY))
MPI_EXECUTABLES = $(addsuffix $(EXECUTABLE_EXT), $(PARALLEL_ARRAY))
EXECUTABLE_PATTERN = %$(EXECUTABLE_EXT)

### Windows/Linux specific setup ###
ifeq ($(OS), Windows_NT)
	MPI_VERSION = 3.1.3

	# additional linking required on Windows
	# Windows environmental variable 'CYGWIN_PATH --> C:\cygwin' required
	MPI_FLAGS += -I "$(CYGWIN_PATH)\usr\include" 
	MPI_FLAGS += -I "$(CYGWIN_PATH)\lib\openmpi"

	EXECUTABLE_EXT = .exe
else
	MPI_VERSION = 3.1.1
endif

### LCC2 specific setup ###
ifeq ($(SGE_CELL), lcc2)
	LOAD_MPI = module load openmpi/$(MPI_VERSION) 2>/dev/null &&

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

	SEQ_EXEC_COMMAND += qsub $(SGE_SEQ_FLAGS) ./$${sequential} $${problem_size}

	MPI_EXEC_COMMAND = echo -e "module load openmpi/$(MPI_VERSION) 2>/dev/null\nmpiexec -n $${num_ranks} ./$${parallel} $${problem_size}" > job.tmp;
	MPI_EXEC_COMMAND += qsub $(SGE_MPI_FLAGS) job.tmp;
	MPI_EXEC_COMMAND += rm job.tmp
else
	SEQ_EXEC_COMMAND = echo "running '$(SEQ_NAME)'";
	SEQ_EXEC_COMMAND += ./$${sequential} $${problem_size} 1> $(SEQ_OUTPUT_FILENAME) 2> $(SEQ_ERROR_FILENAME)

	MPI_EXEC_COMMAND = echo "running '$(MPI_NAME)'";
	MPI_EXEC_COMMAND += mpiexec -n $${num_ranks} --oversubscribe ./$${parallel}$(EXECUTABLE_EXT) $${problem_size} 1> $(MPI_OUTPUT_FILENAME) 2> $(MPI_ERROR_FILENAME)
endif


############## TARGETS ################
.DEFAULT_GOAL := compile_all

### CREATE NON-EXISTING DIRECTORIES ###
$(OUTPUTS_DIR):
	@mkdir -p $@

$(ERRORS_DIR):
	@mkdir -p $@

$(DATA_DIR):
	@mkdir -p $@

$(RESULTS_DIR):
	@mkdir -p $@


### COMPILE ###
%.o: %.$(LANGUAGE) %.h
	@echo "compiling $@"
	@$(CC) $(CC_FLAGS) $(DEFINES) -c $< $(LIBRARIES) -o $@
	@echo

$(SEQ_EXECUTABLES): $(EXECUTABLE_PATTERN): %.$(LANGUAGE) $(DEP_OBJECTS)
	@echo "compiling $@"
	@$(CC) $(CC_FLAGS) $(DEFINES) $^ $(LIBRARIES) -o $@
	@echo

$(MPI_EXECUTABLES): $(EXECUTABLE_PATTERN): %.$(LANGUAGE) $(DEP_OBJECTS)
	@echo "compiling $@"
	@$(LOAD_MPI)$(MPI) $(MPI_FLAGS) $(DEFINES) $^ $(LIBRARIES) -o $@
	@echo

.PHONEY: compile_all
compile_all: $(SEQ_EXECUTABLES) $(MPI_EXECUTABLES)


### RUN ###
.PHONEY: run_seq run_seq_single run_seq_single_clean
run_seq_single_clean: clean run_seq_single

run_seq_single: SEQUENTIAL_ARRAY := $(word 1, $(SEQUENTIAL_ARRAY)) 
run_seq_single: PROBLEM_SIZE_ARRAY := $(word 1, $(PROBLEM_SIZE_ARRAY)) 

run_seq_single: $(word 1, $(SEQ_EXECUTABLES))
run_seq: $(SEQ_EXECUTABLES)

run_seq run_seq_single: $(OUTPUTS_DIR) $(ERRORS_DIR)
	@for sequential in $(SEQUENTIAL_ARRAY); do \
		for problem_size in $(PROBLEM_SIZE_ARRAY); do \
			$(SEQ_EXEC_COMMAND); \
		done \
	done
	@echo

.PHONEY: run_mpi run_mpi_single run_mpi_single_clean
run_mpi_single_clean: clean run_mpi_single

run_mpi_single: PARALLEL_ARRAY := $(word 1, $(PARALLEL_ARRAY))
run_mpi_single: PROBLEM_SIZE_ARRAY := $(word 1, $(PROBLEM_SIZE_ARRAY))
run_mpi_single: NUM_RANKS_ARRAY := $(word 1, $(NUM_RANKS_ARRAY))
run_mpi_single: SLOT_DIRSTRIBUTE_ARRAY := $(word 1, $(SLOT_DISTRIBUTE_ARRAY))

run_mpi_single: $(word 1, $(MPI_EXECUTABLES))
run_mpi: $(MPI_EXECUTABLES)

run_mpi run_mpi_single: $(OUTPUTS_DIR) $(ERRORS_DIR)
	@for parallel in $(PARALLEL_ARRAY); do \
		for run in in $(seq 1 $(NUM_RUNS)); do \
			for num_ranks in $(NUM_RANKS_ARRAY); do \
				for slot_distribute in $(SLOT_DISTRIBUTE_ARRAY); do \
					for problem_size in $(PROBLEM_SIZE_ARRAY); do \
						$(MPI_EXEC_COMMAND); \
					done \
				done \
			done \
		done \
	done
	@echo


.PHONEY: run run_clean
run_clean: clean run

run: run_seq run_mpi
ifeq ($(OS), Windows_NT)
run: check plots
else
	@qstat
endif


### VISUALIZE ###
# TODO --> runtime, speedup, efficiency, ...
.PHONEY: preprocess
preprocess: $(OUTPUTS_DIR) $(DATA_DIR)
	@echo "preprocessing outputs"
	@$(PYTHON) ./$(PREPROCESS) $(OUTPUTS_DIR) $(DATA_DIR)
	@echo

.PHONEY: plots
plots: preprocess $(DATA_DIR) $(RESULTS_DIR)
	@echo "creating some fancy plots..."
	@$(PYTHON) ./$(PLOT) $(DATA_DIR) $(RESULTS_DIR)
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

.PHONEY: check
check: check_errors check_outputs

### CLEANUP ###
# '|| true' --> ignore errors if files don't exist
.PHONEY: clean
clean:
	@rm -r $(ERRORS_DIR) 2>/dev/null  || true
	@rm -r $(OUTPUTS_DIR) 2>/dev/null  || true
	@rm -r $(DATA_DIR) 2>/dev/null  || true
	@rm -r $(RESULTS_DIR) 2>/dev/null  || true
	@echo "your directory has been cleaned errors, outputs, data and results!"
	@echo 

.PHONEY: clean_executables
clean_executables:
	@rm $(DEP_OBJECTS) $(SEQ_EXECUTABLES) $(MPI_EXECUTABLES) 2>/dev/null || true
	@echo "your directory has been cleaned of executables and dependency objects!"
	@echo

.PHONEY: clean_core_dumps
clean_core_dumps:
	@rm *.stackdump 2>/dev/null || true
	@rm core.* 2>/dev/null || true
	@echo "your directory has been cleaned of core dumps!"
	@echo

.PHONEY: super_clean
super_clean: clean clean_executables clean_core_dumps

