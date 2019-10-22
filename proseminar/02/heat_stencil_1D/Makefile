CC=gcc
CC_FLAGS=-O2 -std=c99 -Wall -Wextra -pedantic
MPI=mpicc
Cygwin_Include=-I"C:\cygwin64\usr\include" -I"C:\cygwin64\lib\openmpi"

# sequential
.PHONEY: heat_stencil_1D_seq
heat_stencil_1D_seq: heat_stencil_1D_seq.c
	@$(CC) $(CC_FLAGS) $< -o heat_stencil_1D_seq

.PHONEY: run_seq
run_seq: heat_stencil_1D_seq
	@echo "Sequential:"
	@./heat_stencil_1D_seq
	@echo

# parallel
.PHONEY: heat_stencil_1D_mpi
heat_stencil_1D_mpi: heat_stencil_1D_mpi.c
	@$(MPI) $(CC_FLAGS) $< -o heat_stencil_1D_mpi

.PHONEY: run_mpi
run_mpi: heat_stencil_1D_mpi
	@echo "MPI:"
	@./heat_stencil_1D_mpi
	@echo

# parallel with cygwin
.PHONEY: heat_stencil_1D_mpi_cygwin
heat_stencil_1D_mpi_cygwin: heat_stencil_1D_mpi.c
	@$(MPI) $(CC_FLAGS) $(Cygwin_Include) $< -o heat_stencil_1D_mpi

.PHONEY: run_mpi_cygwin
run_mpi_cygwin: heat_stencil_1D_mpi
	@echo "MPI:"
	@./heat_stencil_1D_mpi
	@echo

# all
.PHONEY: all
all: heat_stencil_1D_seq heat_stencil_1D_mpi

.PHONEY: run
run: run_seq run_mpi

# all with cygwin
.PHONEY: all_cygwin
all_cygwin: heat_stencil_1D_seq heat_stencil_1D_mpi_cygwin

.PHONEY: run_cygwin
run_cygwin: run_seq run_mpi_cygwin


# clean
.PHONEY: clean
clean:
	@rm -f heat_stencil_1D_seq
	@rm -f heat_stencil_1D_mpi
	@rm -f *.o
	
	
	
	