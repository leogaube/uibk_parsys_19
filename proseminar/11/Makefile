CC=gcc
CC_FLAGS=-O3 -std=c99 -g3 -pg -fopenmp
DEFINES=-D _OMP

.PHONEY: all
all: real

SRC=$(wildcard *.c)

real: $(SRC)
	@$(CC) $(CC_FLAGS) $(DEFINES) $(SRC) -lm -o real

.PHONEY: clean
clean:
	@rm real

.PHONEY: run
run: all
	@echo "Sequential:"
	@./real
	@echo
