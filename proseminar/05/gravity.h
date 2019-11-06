#ifndef GRAVITY_H_
#define GRAVITY_H_

#define IDX_2D(x, y, nx) (x)+(y)*(nx)

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

int get_forces(double *forces, Particle_p particles, int N);

#endif /* GRAVITY_H_ */
