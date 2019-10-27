#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

#define RESOLUTION 120

void printTemperature(Vector m, int nx, int ny, int nz);

// -- simulation code ---

int main(int argc, char **argv) {
  clock_t start = clock();

  // 'parsing' optional input parameter = problem size
  int N = 50;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N;
#ifdef VERBOSE
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);
#endif

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N*N*N);

  // set up initial conditions in A
  for (int i = 0; i < N*N*N; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = source_x;
  int source_z = source_z;
  A[IDX_3D(source_x,source_y, source_z,N,N)] = 273 + 60;

#ifdef VERBOSE
  printf("Initial:\n");
  printTemperature(A, N, N, N);
#endif
  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N*N*N);

  // for each time step ..
  for (int t = 0; t < T; t++) {
	// .. we propagate the temperature
	for (int z = 0; z < N; z++) {
		for (int y = 0; y < N; y++) {
			for(int x = 0; x < N; x++){
				// get the current idx
				int i = IDX_2D(x,y,N);

				// center stays constant (the heat is still on)
				if (x == source_x && y == source_y) {
					B[i] = A[i];
					continue;
				}

				// get temperature at current position
				value_t tc = A[i];

				// get temperatures of adjacent cells
				value_t tl = (x != 0) ? A[IDX_3D(x-1,y,z,N,N)] : tc;
				value_t tr = (x != N - 1) ? A[IDX_3D(x+1,y,z,N,N)] : tc;
				value_t tu = (y != 0) ? A[IDX_3D(x,y-1,z,N,N)] : tc;
				value_t td = (y != N - 1) ? A[IDX_3D(x,y+1,z,N,N)] : tc;
				value_t tf = (z != 0) ? A[IDX_3D(x,y,z-1,N,N)] : tc;
				value_t tb = (z != N - 1) ? A[IDX_3D(x,y,z+1,N,N)] : tc;


				// compute new temperature at current position
				B[i] = tc + 0.2 * (tl + tr + tu + td + tf + tb + (-6 * tc));
			}
		}
		// swap matrices (just pointers, not content)
		Vector H = A;
		A = B;
		B = H;
	  }
#ifdef VERBOSE
	printf("time %i:\n",t);
	printTemperature(B, N, N, N);
#endif
  }

  releaseVector(B);
#ifdef VERBOSE
  printf("Final:\n");
  printTemperature(A, N, N, N);
#endif

  // ---------- check ----------
  int success = (is_verified_3D(A, N, N, N, source_x, source_y, source_z, T)==0);
  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  releaseVector(A);

  clock_t end = clock();
  printf("The process took %g seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

void printTemperature(Vector m, int nx, int ny, int nz) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;
  if(W>nx || W>ny || W>nz){
	  W=MIN(MIN(nx,ny),nz);
  }

  // step size in each dimension
  int sWx = nx / W;
  int sWy = ny / W;
  int sWz = nz / W;

  // room
  // actual room
  for(int k=0; k<W; k++){
	  for(int j=0; j<W; j++){
		  // left wall
		  printf("X");
		  for (int i = 0; i < W; i++) {
			// get max temperature in this tile
			value_t max_t = 0;
			for (int z = sWz * k; z < sWz * k + sWz; k++) {
				for (int y = sWy * j; y < sWy * j + sWy; y++) {
					for (int x = sWx * i; x < sWx * i + sWx; x++) {
						max_t = (max_t < m[IDX_3D(x,y,z,nx,ny)]) ? m[IDX_3D(x,y,z,nx,ny)] : max_t;
					}
				}
			}
			value_t temp = max_t;

			// pick the 'color'
			int c = ((temp - min) / (max - min)) * numColors;
			c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

			// print the average temperature
			printf("%c", colors[c]);
		  }
		  // right wall
		  printf("X\n");
	  }
	  printf("\n");
  }
}
