#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

int main (int argc, char** argv) {
	// per default we create 100 million points
	int N = 100000000;
	double start = omp_get_wtime();
    if (argc > 2) {
        N = atoi(argv[1]);
    }
	double pi;
	int count = 0; // number of points in the 1st quadrant of unit circle
	
	// Initialization of random number generator
	#pragma omp parallel for reduction(+:count) 
	for (int i = 0; i < N; i++) {
		unsigned int seed = i^omp_get_thread_num();
		double x = (double) rand_r(&seed) / RAND_MAX;
		seed++;
		double y = (double) rand_r(&seed) / RAND_MAX;
		double z = x * x + y * y;
		if (z <= 1) count++;
	}
	
	pi = (double) count / N * 4;
	double end = omp_get_wtime();
	printf("number of trials = %d, estimate of pi is %g \n", N, pi);
	printf("The process took %g seconds to finish. \n", end - start);

	return EXIT_SUCCESS;
}
