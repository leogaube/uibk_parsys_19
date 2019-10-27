#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

#define RESOLUTION 120

void printTemperature(Vector m, int nx, int ny);

// -- simulation code ---

int main(int argc, char **argv) {
  clock_t start = clock();

  // 'parsing' optional input parameter = problem size
  int N = 50;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = 100;//N * N * 500;
#ifdef VERBOSE
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);
#endif
  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N*N);

  // set up initial conditions in A
  for (int i = 0; i < N*N; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  int source_y = source_x;
  A[IDX_2D(source_x,source_y,N)] = 273 + 60;

#ifdef VERBOSE
  printf("Initial:\n");
  printTemperature(A, N, N);
#endif

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N*N);

  // for each time step ..
  for (int t = 0; t < T; t++) {
	// .. we propagate the temperature
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
//			value_t tl = A[IDX_2D((x-1+N)%N,y,N)];
//			value_t tr = A[IDX_2D((x+1)%N,y,N)];
//			value_t tu = A[IDX_2D(x,(y+1)%N,N)];
//			value_t td = A[IDX_2D(x,(y-1+N)%N,N)];
		      value_t tl = (x != 0) ? A[IDX_2D(x-1,y,N)] : tc;
		      value_t tr = (x != N - 1) ? A[IDX_2D(x+1,y,N)] : tc;
		      value_t tu = (y != 0) ? A[IDX_2D(x,y-1,N)] : tc;
		      value_t td = (y != N - 1) ? A[IDX_2D(x,y+1,N)] : tc;


			// compute new temperature at current position
			B[i] = tc + 0.2 * (tl + tr + tu + td + (-4 * tc));
		}
	}
#ifdef VERBOSE
	printf("time %i:\n",t);
	printTemperature(B, N, N);
#endif
	// swap matrices (just pointers, not content)
	Vector H = A;
	A = B;
	B = H;
  }

  releaseVector(B);
#ifdef VERBOSE
  printf("Final:\n");
  printTemperature(A, N, N);
#endif

  // ---------- check ----------
  int success = (is_verified_2D(A, N, N, source_x, source_y, T)==0);
  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  releaseVector(A);

  clock_t end = clock();
  printf("The process took %g seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

void printTemperature(Vector m, int nx, int ny) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;
  if(W>nx || W>ny){
	  W=MIN(nx,ny);
  }

  // step size in each dimension
  int sWx = nx / W;
  int sWy = ny / W;

  // room
  // actual room
  for(int j=0; j<W; j++){
	  // left wall
	  printf("X");
	  for (int i = 0; i < W; i++) {
		// get max temperature in this tile
		value_t max_t = 0;
		for (int y = sWy * j; y < sWy * j + sWy; y++) {
			for (int x = sWx * i; x < sWx * i + sWx; x++) {
				max_t = (max_t < m[IDX_2D(x,y,nx)]) ? m[IDX_2D(x,y,nx)] : max_t;
			}
		}
		value_t temp = max_t;

		// pick the 'color'
		int c = ((temp - min) / (max - min)) * numColors;
		c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

		// print the average temperature
//		printf("%c", colors[c]);
	    printf("%2.0f\t", m[i+W*j]-273);
	  }
	  // right wall
	  printf("X\n");
  }
}
