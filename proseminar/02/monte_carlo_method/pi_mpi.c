#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main (int argc, char** argv) {
	MPI_Init(&argc, &argv);
	// per default we create 100 million points
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// If N is not a multiple of N than we just dont care because this is negliable in this scenario.
	unsigned long N = 100000000 / size;
	if (argc > 1) {
		sscanf(argv[1], "%lu", &N);
	}
	double x, y, z, pi;
	int i = 0;
	unsigned long count = 0;

	// Initialization of random number generator
	srand(time(NULL));
	for (i = 0; i < N; i++) {
		x = (double) rand() / RAND_MAX;
		y = (double) rand() / RAND_MAX;
		z = x * x + y * y;
		if (z <= 1) count++;
	}
	
	unsigned long pointsInside;
	MPI_Reduce(&count, &pointsInside, 1, MPI_UNSINGED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);


	if (rank == 0) {
		pi = (double) pointsInside / (N * size) * 4;
		printf("number of trials = %lu, estimate of pi is %g \n", N * size, pi);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
