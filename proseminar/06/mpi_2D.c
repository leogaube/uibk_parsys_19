/*
 * mpi_2D.c
 * N-particle simulation with gravity in 2D (parallel)
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <mpi.h>

#include "gravity.h"

// Do not change this, because then the initialization will be fucked up; or change the initialization to ensure the max position is not above or below MAX_POSITION and - MAX_POSITION
#define MAX_POSITION 0.5

int get_forces(double *forces_x, double *forces_y, Particle_p local_particles, int local_rank, Particle_p tmp_particles, int tmp_rank, int M);
int apply_forces(double *forces_x, double *forces_y, Particle_p local_particles, int local_rank, int N, int M);
int init_particles(Particle_p particles, int N, int rank);
int print_particles(Particle_p particles, int N, int room_size);
double let_particles_fly(double position);
int get_com_coords(double* com_coords, Particle_p particles, int N);


int main(int argc, char **argv)
{
	clock_t start = clock();
	int N = 10;
#ifdef VERBOSE
//	int room_size = 10;
#endif
	if (argc == 2){
	    N = atoi(argv[1]);
#ifdef VERBOSE
//		room_size = (N>8) ? 80 : N*10;
#endif
	}
	if (argc == 3) {
		N = atoi(argv[1]);
#ifdef VERBOSE
//		room_size = atoi(argv[2]);
#endif
	}
	int T = 10;

    // MPI setup
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (N % numProcs != 0){
        printf("This Problem cannot be split up evenly among MPI ranks! (N mod numProcs != 0)");
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    int M = N / numProcs;

    MPI_Datatype twoDVector;
    int blocklength[1] = {2};
    MPI_Aint displacement[1] = {offsetof(TwoDVector,x)};
    MPI_Datatype datatype[1] = {MPI_DOUBLE};
    MPI_Type_create_struct(1, blocklength, displacement, datatype, &twoDVector);
    MPI_Type_commit(&twoDVector);

    MPI_Datatype particles_type;
    int blocklengths[2] = {1,2};
    MPI_Aint displacements[2] = {offsetof(Particle,mass), offsetof(Particle,velocity)};
    MPI_Datatype datatypes[2] = {MPI_DOUBLE, twoDVector};
    MPI_Type_create_struct(2, blocklengths, displacements, datatypes, &particles_type);
    MPI_Type_commit(&particles_type);

    // init particles with random values
	Particle_p local_particles = malloc(M*sizeof(Particle));
	init_particles(local_particles, M, rank);
	
	// this particle buffer will be loaded with the particles from the other ranks
	Particle_p tmp_particles = malloc(M*sizeof(Particle));

	// where are the tmp_particles from (in the first step there has not been any communication)
	int tmp_particles_rank = rank;

	// only upper triangular matrix (without diagonal) is needed
	double *forces_x = malloc(((int)(N*(N-1)/2))*sizeof(double));
	double *forces_y = malloc(((int)(N*(N-1)/2))*sizeof(double));
	for(int t=0; t<T; t++){
		for(int i=0; i<numProcs; i++){
			// calculate the forces between the local_particles and the tmp_particles
			get_forces(forces_x, forces_y, local_particles, rank, tmp_particles, tmp_particles_rank, M);

			// send/recveive the local particles around
			if(i<numProcs-1){
				tmp_particles_rank = (rank+i+1)%numProcs;
				if(rank%2==0){
					MPI_Ssend(local_particles, M, particles_type, (rank-(i+1))%numProcs, 0, MPI_COMM_WORLD);
					MPI_Recv(tmp_particles, M, particles_type, tmp_particles_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				} else {
					MPI_Recv(tmp_particles, M, particles_type, tmp_particles_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Ssend(local_particles, M, particles_type, (rank-(i+1))%numProcs, 0, MPI_COMM_WORLD);
				}
			}

			//TODO possible optimization send/receive forces
		}
		apply_forces(forces_x, forces_y, local_particles, rank, N, M);
	}

	free(forces_x);
	free(forces_y);
	free(local_particles);
	free(tmp_particles);
    if(rank == 0) {
#ifdef VERBOSE
//    	TODO gather all particles for plotting
//		print_particles(particles, N, room_size);
#endif
        clock_t end = clock();
	    printf("The process took %f seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);
    }
	
  	MPI_Finalize();

	return EXIT_SUCCESS;
}

/**
 * calculates the center of mass for verification purposes
 */
int get_com_coords(double* com_coords, Particle_p particles, int N){
	double M = 0;

	for(int i=0; i<N; i++){
		Particle pi = particles[i];
		com_coords[0] += pi.mass*pi.position.x;
		com_coords[1] += pi.mass*pi.position.y;
		M += pi.mass;
	}
	com_coords[0] /= M;
	com_coords[1] /= M;

	return EXIT_SUCCESS;
}

/**
 * calculates the forces in x and y direction between the N particles
 * only the upper triangular matrix (without diagonal) is calculated (i>j)
 */
int get_forces(double *forces_x, double *forces_y, Particle_p local_particles, int local_rank, Particle_p tmp_particles, int tmp_rank, int M){

	if(local_rank == tmp_rank) { //first call -> only calculate between own particles
		for(int i=1; i<M; i++){
			Particle pi = local_particles[i];
			double mi = pi.mass;
			double xi = pi.position.x;
			double yi = pi.position.y;

			for(int j=0; j<i; j++){
				Particle pj = local_particles[j];
				double mj = pj.mass;
				double dx = xi - pj.position.x;
				double dy = yi - pj.position.y;

				// calculate the force in the given directions
				double tmp = - mi*mj/pow(dx*dx+dy*dy,1.5);
				forces_x[IDX_FORCES_MPI(i, j, local_rank, local_rank, M)] = tmp*dx;
				forces_y[IDX_FORCES_MPI(i, j, local_rank, local_rank, M)] = tmp*dy;
			}
		}
	} else { // calculate all forces between the local and external particles
		for(int i=0; i<M; i++){
			Particle pi = (local_rank>tmp_rank) ? local_particles[i] : tmp_particles[i];
			double mi = pi.mass;
			double xi = pi.position.x;
			double yi = pi.position.y;

			for(int j=0; j<M; j++){
				Particle pj = (local_rank>tmp_rank) ? tmp_particles[j] : local_particles[j];
				double mj = pj.mass;
				double dx = xi - pj.position.x;
				double dy = yi - pj.position.y;

				// calculate the force in the given directions
				double tmp = - mi*mj/pow(dx*dx+dy*dy,1.5);
				int idx = (local_rank>tmp_rank) ? IDX_FORCES_MPI(i, j, local_rank, tmp_rank, M) : IDX_FORCES_MPI(i, j, tmp_rank, local_rank, M);
				forces_x[idx] = tmp*dx;
				forces_y[idx] = tmp*dy;
			}
		}
	}


	return EXIT_SUCCESS;
}


/**
 * the total force on each particle is the superposition of all forces
 * calculate this sum and apply it to the position and velocity of the particles
 */
int apply_forces(double *forces_x, double *forces_y, Particle_p local_particles, int local_rank, int N, int M){
	for(int i=0; i<M; i++){
		// get total forces
		double force_x = 0;
		double force_y = 0;

		for(int j=0; j<N; j++){
			force_x += (j < i+M*local_rank) ? forces_x[IDX_FORCES(i+M*local_rank,j)] : -forces_x[IDX_FORCES(j,i+M*local_rank)];
			force_y += (j < i+M*local_rank) ? forces_y[IDX_FORCES(i+M*local_rank,j)] : -forces_y[IDX_FORCES(j,i+M*local_rank)];
		}
		double m = local_particles[i].mass;
		local_particles[i].velocity.x += force_x/m;
		local_particles[i].velocity.y += force_y/m;
		local_particles[i].position.x += local_particles[i].velocity.x;
		local_particles[i].position.y += local_particles[i].velocity.y;
		if (local_particles[i].position.x > MAX_POSITION || local_particles[i].position.x < -MAX_POSITION) {
			local_particles[i].position.x = let_particles_fly(local_particles[i].position.x);
		}
		if (local_particles[i].position.y > MAX_POSITION || local_particles[i].position.y < -MAX_POSITION) {
			local_particles[i].position.y = let_particles_fly(local_particles[i].position.y);
		}
	}

	return EXIT_SUCCESS;
}

/**
 * initializes the particles with random mass and random position
 */
int init_particles(Particle_p particles, int N, int rank) {
	srand(time(NULL) ^ rank);

	for (int i = 0; i < N; i++) {
		particles[i].mass = (rand() / (double) RAND_MAX * 0.9 + 0.1) * 1e-2/N;
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
			// c gives how many mean particle masses are within a cell rounded up
			int c = (int) ceil((room_printed[i][j] / 0.005 * N ));
			c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);
			printf("%c", colors[c]);
		}
		printf("X\n");
	}

	return EXIT_SUCCESS;
}

double let_particles_fly(double position) {
	if (position < -MAX_POSITION) {
		while (position < -MAX_POSITION) {
			position += MAX_POSITION*2;
		}
		return position;
	} else if (position > MAX_POSITION) {
		while (position > MAX_POSITION) {
			position -= MAX_POSITION*2;
		}
		return position;
	}
	return position;
}
