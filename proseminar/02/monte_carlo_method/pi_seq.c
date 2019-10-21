#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

int main (int argc, char** argv) {
	clock_t start = clock();
	// per default we create 100 million points
	int N = 100000000;
	if (argc > 1) {
		N = atoi(argv[1]);
	}
	double x, y, z, pi;
	int i, count = 0; // number of points in the 1st quadrant of unit circle
	
	// Initialization of random number generator
	srand(time(NULL));
	for (i = 0; i < N; i++) {
		x = (double) rand() / RAND_MAX;
		y = (double) rand() / RAND_MAX;
		z = x * x + y * y;
		if (z <= 1) count++;
	}
	
	pi = (double) count / N * 4;
	clock_t end = clock();
	printf("number of trials = %d, estimate of pi is %g \n", N, pi);
	printf("The process took %g seconds to finish. \n",((double) (end - start))/CLOCKS_PER_SEC);

	return EXIT_SUCCESS;
}
