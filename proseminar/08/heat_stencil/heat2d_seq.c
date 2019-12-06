#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

// -- simulation code ---

int main(int argc, char **argv) {
  clock_t start = clock();

  // 'parsing' optional input parameter = problem size
  int Nx = 10;
  int Ny = 10;
  if (argc == 2) {
    Nx = Ny = atoi(argv[1]);
  } else if (argc == 3){
	  Nx = atoi(argv[1]);
	  Ny = atoi(argv[2]);
  }
  int T = MAX(Nx, Ny) * 100;
#ifdef VERBOSE
  printf("Computing heat-distribution for room size Nx=%d, Ny=%d for T=%d timesteps\n", Nx, Ny, T);
#endif
  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx*Ny);

  // set up initial conditions in A
  for (long i = 0; i < Nx*Ny; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  A[IDX_2D(source_x,source_y,Nx)] = 273 + 60;

#ifdef VERBOSE
  printf("Initial:\n");
  printTemperature(A, Nx, Ny, 1);
#endif

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx*Ny);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (int y = 0; y < Ny; y++)
    {
      for (int x = 0; x < Nx; x++)
      {
        // get the current idx
        long i = IDX_2D(x, y, Nx);

        // center stays constant (the heat is still on)
        if (x == source_x && y == source_y)
        {
          B[i] = A[i];
          continue;
        }

        // get temperature at current position
        value_t tc = A[i];

        // get temperatures of adjacent cells
        value_t tl = (x != 0) ? A[IDX_2D(x - 1, y, Nx)] : tc;
        value_t tr = (x != Nx - 1) ? A[IDX_2D(x + 1, y, Nx)] : tc;
        value_t tu = (y != 0) ? A[IDX_2D(x, y - 1, Nx)] : tc;
        value_t td = (y != Ny - 1) ? A[IDX_2D(x, y + 1, Nx)] : tc;

        // compute new temperature at current position
        B[i] = tc + 0.4/4 * (tl + tr + tu + td + (-4 * tc));
      }
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

#ifdef VERBOSE
    if (!(t % 1000))
    {
      printf("Step t=%d:\n", t);
      printTemperature(A, Nx, Ny, 1);
      printf("\n");
    }
#endif
  }

  releaseVector(B);
#ifdef VERBOSE
  printf("Final:\n");
  printTemperature(A, Nx, Ny, 1);
#endif

  // ---------- check ----------
  double residual = is_verified_2D(A, Nx, Ny, source_x, source_y, T);
  printf("The maximal deviation from the 1D theory is %fK.", residual);

  // ---------- cleanup ----------

  releaseVector(A);

  clock_t end = clock();
  printf("The process took %f seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);

  // done
  return EXIT_SUCCESS;
}
