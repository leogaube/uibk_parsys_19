/*
 * seq_2D.c
 * N-particle simulation with gravity in 2D (sequential)
 */
#include <stdio.h>
#include <stdlib.h>

#include "gravity.h"

int main(int argc, char **argv)
{
	int N = 10;
	if (argc == 2){
	    N = atoi(argv[1]);
	}
	int T = N*10;

	// TODO init particles with random values
	Particle_p particles = malloc(N*sizeof(Particle));

	// only upper triangular matrix (without diagonal) is needed
	double *forces = malloc(((int)(N*(N-1)/2))*sizeof(double));
	for(int t=0; t<T; t++){
		get_forces(forces, particles, N);

		// TODO apply forces

		//TODO plot results
	}

	// TODO verification

	free(forces);
	free(particles);
}
