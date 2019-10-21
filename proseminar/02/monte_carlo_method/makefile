CC=gcc
CC_FLAGS=-O2 -std=c99 -Wall
MPI=mpicc
Cygwin_Include=-I"C:\cygwin64\usr\include" -I"C:\cygwin64\lib\openmpi"

# sequential
.PHONEY: seq
seq: pi_seq

pi_seq: pi_seq.c
	@$(CC) $(CC_FLAGS) $< -o pi_seq
	
.PHONEY: run_seq
run_seq: seq
	@echo "Sequential:"
	@./pi_seq
	@echo
	

# parallel
.PHONEY: mpi
mpi: pi_mpi

pi_mpi: pi_mpi.c
	@$(MPI) $(CC_FLAGS) $< -o pi_mpi

.PHONEY: run_mpi
run_mpi: mpi
	@echo "MPI:"
	@./pi_mpi
	@echo


# parallel with cygwin
.PHONEY: mpi_cygwin
mpi_cygwin: pi_mpi

pi_mpi_cygwin: pi_mpi.c
	@$(MPI) $(CC_FLAGS) $(Cygwin_Include) $< -o pi_mpi

.PHONEY: run_mpi_cygwin
run_mpi_cygwin: mpi_cygwin
	@echo "MPI:"
	@./pi_mpi
	@echo
		
	
# all
.PHONEY: all
all: seq mpi

.PHONEY: run
run: run_seq run_mpi


# all with cygwin
.PHONEY: all_cygwin
all_cygwin: seq mpi_cygwin

.PHONEY: run_cygwin
run_cygwin: run_seq run_mpi_cygwin


# clean
.PHONEY: clean
clean:
	@rm -f pi_seq
	@rm -f pi_mpi
	@rm -f *.o
