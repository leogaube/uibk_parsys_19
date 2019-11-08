/*
 * seq_2D.c
 * N-particle simulation with gravity in 2D (sequential)
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

#include "gravity.h"

// Do not change this, because then the initialization will be fucked up ;)
#define MAX_POSITION 0.5

int get_forces(double *forces_x, double *forces_y, Particle_p particles, int N);
int apply_forces(double *forces_x, double *forces_y, Particle_p particles, int N);
int init_particles(Particle_p particles, int N);
int print_particles(Particle_p Particles, int N, int room_size);


int main(int argc, char **argv)
{
	int N = 10;
	int room_size = 10;
	if (argc == 2){
	    N = atoi(argv[1]);
		room_size = N*10;
	}
	if (argc == 3) {
		N = atoi(argv[1]);
		room_size = atoi(argv[2]);
	}
	int T = N*10;

	// init particles with random values
	Particle_p particles = malloc(N*sizeof(Particle));
	init_particles(particles, N);
	
	#ifdef VERBOSE
	printf("Init:\n");
	print_particles(particles, N, room_size);
	#endif

	// only upper triangular matrix (without diagonal) is needed
	double *forces_x = malloc(((int)(N*(N-1)/2))*sizeof(double));
	double *forces_y = malloc(((int)(N*(N-1)/2))*sizeof(double));
	for(int t=0; t<T; t++){
		get_forces(forces_x, forces_y, particles, N);

		apply_forces(forces_x, forces_y, particles, N);

		// plot results
		#ifdef VERBOSE
		if (t % 1 == 0) {
			printf("time t: %d\n", t);
			print_particles(particles, N, room_size);
		}
		#endif
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
			if(i==j){ // there is no force from the particle itself
				continue;
			}
			force_x += (j<i) ? forces_x[IDX_FORCES(i,j)] : -forces_x[IDX_FORCES(j,i)];
			force_y += (j<i) ? forces_y[IDX_FORCES(i,j)] : -forces_y[IDX_FORCES(j,i)];
		}
		double m = particles[i].mass;
		particles[i].velocity.x += force_x/m;
		particles[i].velocity.y += force_y/m;
		particles[i].position.x += particles[i].velocity.x;
		particles[i].position.y += particles[i].velocity.y;
		if (particles[i].position.x > MAX_POSITION) {
			particles[i].position.x = -0.5 + fmod(particles[i].position.x, MAX_POSITION);
		}
		if (particles[i].position.x < -MAX_POSITION) {
			particles[i].position.x = 0.5 + fmod(particles[i].position.x, MAX_POSITION);
		}
		if (particles[i].position.y > MAX_POSITION) {

			particles[i].position.y = -0.5 + fmod(particles[i].position.y, MAX_POSITION);
		}
		if (particles[i].position.y < -MAX_POSITION) {
			particles[i].position.y = 0.5 + fmod(particles[i].position.y, MAX_POSITION);
		}
	}

	return EXIT_SUCCESS;
}

/**
 * initializes the particles with random mass and random position
 */
int init_particles(Particle_p particles, int N) {
	srand(time(NULL));

	for (int i = 0; i < N; i++) {
		// TODO what random values for mass should be generated?
		particles[i].mass = /*(rand() / (double) RAND_MAX * 0.9 + 0.1) **/ 1e-2/N;
		particles[i].position.x = (rand() / (double) ((unsigned)RAND_MAX + 1) - MAX_POSITION);
		particles[i].position.y = (rand() / (double) ((unsigned)RAND_MAX + 1) - MAX_POSITION);
		particles[i].velocity.x = 0;
		particles[i].velocity.y = 0;
	}

	return EXIT_SUCCESS;
}

/**
 * function for checking of a particle is inside a specified index in relation to the room size
 */
double get_mass_in_index(Particle_p particles, int N, int i, int j, int room_size) {
	double m = 0.0;
	for (int k = 0; k < N; k++) {
		if ((int) ((particles[k].position.x + MAX_POSITION) * room_size) == j && (int) ((particles[k].position.y + MAX_POSITION) * room_size) == i) {
			m += particles[k].mass;
		}
	}
	return m;
}

/**
 * function for printing the room where the particles move
 * the movement of the particles is in relation to the room size.
 */
int print_particles(Particle_p particles, int N, int room_size) {
	double room_printed[room_size][room_size];
	const char *colors = " .-:=+*^X#%@";
	const int numColors = 12;

	for (int i = 0; i < room_size; i++) {
		for (int j = 0; j < room_size; j++) {
			room_printed[i][j] = get_mass_in_index(particles, N, i, j, room_size);
		}
	}

	for (int i = 0; i < room_size; i++) {
		printf("X");
		for (int j = 0; j < room_size; j++) {
			int c = (room_printed[i][j] / 0.01) * numColors;
			c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);
			printf("%c", colors[c]);
			//printf("%2.4f ", room_printed[i][j]);//colors[c]);
		}
		printf("X\n");
	}

	return EXIT_SUCCESS;
}
