CC=gcc
CC_FLAGS=-O3 -g3 -lm -pg
OMP_FLAGS=$(CC_FLAGS) -fopenmp

.PHONEY: all
all: real real_omp

SRC=$(wildcard *.c)

real: $(SRC)
	@$(CC) $(CC_FLAGS) $(SRC) -o real
 
real_omp: $(SRC)
	@$(CC) $(OMP_FLAGS) $(SRC) -o real_omp

.PHONEY: clean
clean:
	@rm real || true
	@rm real_omp || true

.PHONEY: run
run: all
	@echo "OMP:"
	@./real_omp
	@echo
	@echo "Sequential:"
	@./real
	@echo
