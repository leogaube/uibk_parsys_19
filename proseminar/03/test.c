/*
 * test.c
 *
 *  Created on: 27.10.2019
 *      Author: marc_bussjaeger
 */
#include <stdio.h>
#include <stdlib.h>

#define IDX_3D(x,y,z,nx,ny) (x+y*nx+z*nx*ny)
#define IDX_2D(x,y,nx) IDX_3D(x,y,0,nx,0)

typedef double value_t;

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);


int main(int argc, char **argv) {
  int N = 2000;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int nx = 2;
  int ny = 3;
  int nz = 1;
  int i = 0;
  Vector A = createVector(nx*ny*nz);
  for(int z=0; z<nz; z++){
	  for(int y=0; y<ny; y++){
		  for(int x=0; x<nx; x++){
			  A[IDX_2D(x,y,nx)] = i;
			  i++;
		  }
	  }
  }
  for(int z=0; z<nz; z++){
	  for(int y=0; y<ny; y++){
		  for(int x=0; x<nx; x++){
			  printf("%f ", A[IDX_2D(x,y,nx)]);
		  }
		  printf("-----\n");
	  }
	  printf("\n");
  }
  printf("\n");
  for(int j=0; j<nx*ny*nz; j++){
	  printf("%f ", A[j]);
  }
  printf("\n");
  releaseVector(A);

  return EXIT_SUCCESS;
}



Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

