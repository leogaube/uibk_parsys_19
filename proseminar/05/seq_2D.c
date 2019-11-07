/*
 * seq_2D.c
 * N-particle simulation with gravity in 2D (sequential)
 */
#include <stdio.h>
#include <stdlib.h>

#include "gravity.h"

int get_forces(double *forces_x, double *forces_y, Particle_p particles, int N);
int apply_forces(double *forces_x, double *forces_y, Particle_p particles, int N);


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
	double *forces_x = malloc(((int)(N*(N-1)/2))*sizeof(double));
	double *forces_y = malloc(((int)(N*(N-1)/2))*sizeof(double));
	for(int t=0; t<T; t++){
		get_forces(forces_x, forces_y, particles, N);

		apply_forces(forces_x, forces_y, particles, N);

		//TODO plot results
	}

	// TODO verification

	free(forces_x);
	free(forces_y);
	free(particles);
}


/**
 * calculates the forces in x and y direction between the N particles
 * only the upper triangular matrix (without diagonal) is calculated (i>j)
 */
int get_forces(double *forces_x, double *forces_y, Particle_p particles, int N){

	for(int i=1; i<N; i++){
		Particle pi = particles[i];
		double mi = pi.mass;
		double xi = pi.position.x;
		double yi = pi.position.y;

		for(int j=0; j<i; j++){
			Particle pj = particles[j];
			double mj = pj.mass;
			double dx = xi - pj.position.x;
			double dy = yi - pj.position.y;

			// calculate the force in the given directions
			double tmp = - mi*mj/pow(dx*dx+dy*dy,1.5);
			forces_x[IDX_FORCES(i, j)] = tmp*dx;
			forces_y[IDX_FORCES(i, j)] = tmp*dy;
		}
	}

	return EXIT_SUCCESS;
}


/**
 * the total force on each particle is the superposition of all forces
 * calculate this sum and apply it to the position and velocity of the particles
 */
int apply_forces(double *forces_x, double *forces_y, Particle_p particles, int N){
	for(int i=0; i<N; i++){
		// get total forces
		double force_x = 0;
		double force_y = 0;
		for(int j=0; j<N; j++){
			force_x += (j<i) ? forces_x[IDX_FORCES(i,j)] : -forces_x[IDX_FORCES(j,i)];
			force_y += (j<i) ? forces_y[IDX_FORCES(i,j)] : -forces_y[IDX_FORCES(j,i)];
		}
		double m = particles[i].mass;
		particles[i].velocity.x += force_x/m;
		particles[i].velocity.y += force_y/m;
		particles[i].position.x += particles[i].velocity.x;
		particles[i].position.y += particles[i].velocity.y;
	}

	return EXIT_SUCCESS;
}



