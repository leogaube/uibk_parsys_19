#ifndef HEAT_STENCIL_H_
#define HEAT_STENCIL_H_

typedef float value_t;

#define RESOLUTION 25

#define IDX_3D(x,y,z,nx,ny) ((x)+(y)*(nx)+(z)*(nx)*(ny))
#define IDX_2D(x,y,nx) IDX_3D(x,y,0,nx,0)

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int nx, int ny, int nz);

// -- verification --
double is_verified_2D(Vector result_2D, int nx, int ny, int source_x, int source_y, int T);
double is_verified_3D(Vector result_3D, int nx, int ny, int nz, int source_x, int source_y, int source_z, int T);



#endif /* HEAT_STENCIL_H_ */
