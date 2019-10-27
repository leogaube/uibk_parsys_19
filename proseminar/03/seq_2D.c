#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

#define RESOLUTION 120


// -- simulation code ---

int main(int argc, char **argv) {
  clock_t start = clock();

  // 'parsing' optional input parameter = problem size
  int N = 50;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N * 500;
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

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
			value_t tl = (x != 0) ? A[IDX_2D(x-1,y,N)] : tc;
			value_t tr = (x != N - 1) ? A[IDX_2D(x+1,y,N)] : tc;
			value_t tu = (y != 0) ? A[IDX_2D(x,y-1,N)] : tc;
			value_t td = (y != N - 1) ? A[IDX_2D(x,y+1,N)] : tc;

			// compute new temperature at current position
			B[i] = tc + 0.2 * (tl + tr + tu + td + (-4 * tc));
		}

		// swap matrices (just pointers, not content)
		Vector H = A;
		A = B;
		B = H;
	}
  }

  releaseVector(B);

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

