#ifndef GRAVITY_H_
#define GRAVITY_H_

// using the upper triangular matrix the index of i has to be larger than j
#define IDX_FORCES(i, j) (j)+(i)*(i-1)/2

typedef struct{
	double x;
	double y;
} TwoDVector;

typedef struct{
	double mass;
	TwoDVector velocity;
	TwoDVector position;
} Particle;

typedef Particle *Particle_p;

#endif /* GRAVITY_H_ */
