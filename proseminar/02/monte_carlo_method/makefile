CC=gcc
CC_FLAGS=-O2 -std=c99 -Wall
MPI=mpicc

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
	
	
# all
.PHONEY: all
all: seq mpi

.PHONEY: run
run: run_seq run_mpi


# clean
.PHONEY: clean
clean:
	@rm pi_seq
	@rm pi_mpi
	@rm *.o
