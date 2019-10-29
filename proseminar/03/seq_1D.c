#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

// -- simulation code ---

int main(int argc, char **argv) {
  clock_t start = clock();

  // 'parsing' optional input parameter = problem size
  int N = 2000;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N * 500;
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N);

  // set up initial conditions in A
  for (int i = 0; i < N; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  A[source_x] = 273 + 60;

  printf("Initial:\t");
  printTemperature(A, N, 1, 1);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (long long i = 0; i < N; i++) {
      // center stays constant (the heat is still on)
      if (i == source_x) {
        B[i] = A[i];
        continue;
      }

      // get temperature at current position
      value_t tc = A[i];

      // get temperatures of adjacent cells
      value_t tl = (i != 0) ? A[i - 1] : tc;
      value_t tr = (i != N - 1) ? A[i + 1] : tc;

      // compute new temperature at current position
      B[i] = tc + 0.16666 * (tl + tr + (-2 * tc));
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 1000)) {
      printf("Step t=%d:\t", t);
      printTemperature(A, N, 1, 1);
      printf("\n");
    }
  }

  releaseVector(B);

  // ---------- check ----------

  printf("Final:\t\t");
  printTemperature(A, N, 1, 1);
  printf("\n");

  int success = 1;
  for (long long i = 0; i < N; i++) {
    value_t temp = A[i];
    if (273 <= temp && temp <= 273 + 60)
      continue;
    success = 0;
    break;
  }

  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  releaseVector(A);

  clock_t end = clock();
  printf("The process took %f seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}