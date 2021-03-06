/*
 * omp_2D.c
 * N-particle simulation with gravity in 2D (with OpenMP)
 */
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "gravity.h"

// Do not change this, because then the initialization will be fucked up; or change the initialization to ensure the max position is not above or below MAX_POSITION and - MAX_POSITION
#define MAX_POSITION 0.5

int get_forces(double *forces_x, double *forces_y, Particle_p particles, int N);
int apply_forces(double *forces_x, double *forces_y, Particle_p particles, int N);
int init_particles(Particle_p particles, int N);
int print_particles(Particle_p particles, int N, int room_size);
double let_particles_fly(double position);
int get_com_coords(double *com_coords, Particle_p particles, int N);

int main(int argc, char **argv) {
    omp_set_num_threads(6);

    double start = omp_get_wtime();
    int N = 3000;
    int room_size = 100;
    if (argc == 2) {
        N = atoi(argv[1]);
        room_size = (N > 8) ? 80 : N * 10;
    }
    if (argc == 3) {
        N = atoi(argv[1]);
        room_size = atoi(argv[2]);
    }
    int T = 100;

    // init particles with random values
    Particle_p particles = malloc(N * sizeof(Particle));
    init_particles(particles, N);

    // get the center of mass coordinates
    double com_coords[2];
    get_com_coords(com_coords, particles, N);

    // only upper triangular matrix (without diagonal) is needed
    double *forces_x = malloc(((int)(N * (N - 1) / 2)) * sizeof(double));
    double *forces_y = malloc(((int)(N * (N - 1) / 2)) * sizeof(double));

    for (int t = 0; t < T; t++) {
        get_forces(forces_x, forces_y, particles, N);
        apply_forces(forces_x, forces_y, particles, N);
    }
    double end = omp_get_wtime();

#ifdef VERBOSE
    print_particles(particles, N, room_size);
#endif

    printf("The process took %f seconds to finish. \n", (end - start));

    // verification
    double com_coords_T[2];
    get_com_coords(com_coords_T, particles, N);
    printf("The COM moved by (%f,%f) units\n", com_coords_T[0] - com_coords[0], com_coords_T[1] - com_coords[1]);

    free(forces_x);
    free(forces_y);
    free(particles);

    return EXIT_SUCCESS;
}

/**
 * calculates the center of mass for verification purposes
 */
int get_com_coords(double *com_coords, Particle_p particles, int N) {
    double M = 0;

    for (int i = 0; i < N; i++) {
        Particle pi = particles[i];
        com_coords[0] += pi.mass * pi.position.x;
        com_coords[1] += pi.mass * pi.position.y;
        M += pi.mass;
    }
    com_coords[0] /= M;
    com_coords[1] /= M;

    return EXIT_SUCCESS;
}

int get_i_from_N(int N){
	return (int) floor(-1/2+sqrt(1/4+2*N))+1;
}

int get_j_from_N(int N, int i){
	return N-i*(i-1)/2;
}

/**
 * calculates the forces in x and y direction between the N particles
 * only the upper triangular matrix (without diagonal) is calculated (i>j)
 */
int get_forces(double *forces_x, double *forces_y, Particle_p particles, int N) {
	Particle pi;
	int i_tmp = 0;
#pragma omp parallel for
	for (int k = 0; k < N*(N-1)/2; k++) {
    	int i = get_i_from_N(k);
        //printf("%d\n", omp_get_num_threads());
		double mi;
		double xi;
		double yi;
        if(i!=i_tmp){
        	pi = particles[i];
			mi = pi.mass;
			xi = pi.position.x;
			yi = pi.position.y;
			i_tmp=i;
        }

        int j = get_j_from_N(k, i);
		Particle pj = particles[j];
		double mj = pj.mass;
		double dx = xi - pj.position.x;
		double dy = yi - pj.position.y;

		// calculate the force in the given directions
		double tmp = -mi * mj / pow(dx * dx + dy * dy, 1.5);
		forces_x[IDX_FORCES(i, j)] = tmp * dx;
		forces_y[IDX_FORCES(i, j)] = tmp * dy;

//        if(IDX_FORCES(i,j)>=N*(N-1)/2){
//        	printf("idx too large(%i>=%i) i=%i, j=%i, k=%i\n",IDX_FORCES(i,j),N*(N-1)/2,i,j,k);
//        }
//        if(j>=i){
//        	printf("i too small i=%i, j=%i, k=%i\n",i,j,k);
//        }
//        if(IDX_FORCES(i,j)!=k){
//        	printf("i,j does not give k but k+%i\n",IDX_FORCES(i,j)+k);
//        }

    }

    return EXIT_SUCCESS;
}

/**
 * the total force on each particle is the superposition of all forces
 * calculate this sum and apply it to the position and velocity of the particles
 */
int apply_forces(double *forces_x, double *forces_y, Particle_p particles, int N) {
#pragma omp parallel for
    for (int i = 0; i < N; i++) {
        // get total forces
        double force_x = 0;
        double force_y = 0;
        for (int j = 0; j < N; j++) {
            if (i == j) {  // there is no force from the particle itself
                continue;
            }
            force_x += (j < i) ? forces_x[IDX_FORCES(i, j)] : -forces_x[IDX_FORCES(j, i)];
            force_y += (j < i) ? forces_y[IDX_FORCES(i, j)] : -forces_y[IDX_FORCES(j, i)];
        }
        double m = particles[i].mass;
        particles[i].velocity.x += force_x / m;
        particles[i].velocity.y += force_y / m;
        particles[i].position.x += particles[i].velocity.x;
        particles[i].position.y += particles[i].velocity.y;
        if (particles[i].position.x > MAX_POSITION || particles[i].position.x < -MAX_POSITION) {
            particles[i].position.x = let_particles_fly(particles[i].position.x);
        }
        if (particles[i].position.y > MAX_POSITION || particles[i].position.y < -MAX_POSITION) {
            particles[i].position.y = let_particles_fly(particles[i].position.y);
        }
    }

    return EXIT_SUCCESS;
}

/**
 * initializes the particles with random mass and random position
 */
int init_particles(Particle_p particles, int N) {
    srand(1234);

    for (int i = 0; i < N; i++) {
        particles[i].mass = (rand() / (double)RAND_MAX * 0.9 + 0.1) * 1e-2 / N;
        particles[i].position.x = (rand() / (double)((unsigned)RAND_MAX + 1) - MAX_POSITION);
        particles[i].position.y = (rand() / (double)((unsigned)RAND_MAX + 1) - MAX_POSITION);
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
        if ((int)((particles[k].position.x + MAX_POSITION) * room_size) == j && (int)((particles[k].position.y + MAX_POSITION) * room_size) == i) {
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
            int c = (int)ceil((room_printed[i][j] / 0.005 * N));
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
            position += MAX_POSITION * 2;
        }
        return position;
    } else if (position > MAX_POSITION) {
        while (position > MAX_POSITION) {
            position -= MAX_POSITION * 2;
        }
        return position;
    }
    return position;
}
