CC=gcc
CC_FLAGS=-O2 -std=c99 -Wall
MPI=mpicc
MPI_FLAGS=-O2 -std=c99 -Wall

#on windows add include directories 
ifeq ($(OS), Windows_NT)
	NUM_RANKS=4
	MPI_FLAGS+=-I"C:\cygwin64\usr\include" -I"C:\cygwin64\lib\openmpi"
endif

# sequential
.PHONEY: heat_stencil_1D_seq
heat_stencil_1D_seq: heat_stencil_1D_seq.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_1D_seq

.PHONEY: run_seq
run_seq: heat_stencil_1D_seq
	@echo "Sequential:"
ifeq ($(OS), Windows_NT)
	@./heat_stencil_1D_seq
else
	@qsub job_seq.script
endif
	@echo

# parallel
.PHONEY: heat_stencil_1D_mpi Makefile
heat_stencil_1D_mpi: heat_stencil_1D_mpi.c
	@$(MPI) $(MPI_FLAGS) $< -o heat_stencil_1D_mpi

.PHONEY: run_mpi
run_mpi: heat_stencil_1D_mpi
	@echo "MPI:"
ifeq ($(OS), Windows_NT)
	@mpiexec -n $(NUM_RANKS) heat_stencil_1D_mpi
else
	@qsub job_fillup.script
	@qsub job_2perhost.script
endif
	@echo

# all
.PHONEY: all
all: heat_stencil_1D_seq heat_stencil_1D_mpi

.PHONEY: run
run: run_seq run_mpi

# clean
.PHONEY: clean
clean:
	@rm heat_stencil_1D_seq 2> /dev/null || true
	@rm heat_stencil_1D_mpi 2> /dev/null || true
	@rm output* 2> /dev/null || true
