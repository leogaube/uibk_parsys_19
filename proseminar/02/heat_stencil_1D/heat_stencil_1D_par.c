#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 120

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int N);

// -- simulation code ---

int main(int argc, char **argv) {
// 'parsing' optional input parameter = problem size
  int N = 2000;
  if (argc > 1)
  {
    N = atoi(argv[1]);
  }
  int T = N * 500;
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

  // MPI setup
  int myrank, numProcs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (N % numProcs != 0){
    printf("This Problem cannot be split up evenly among MPI ranks! (N mod numProcs != 0)");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  int elemsPerRank = N / numProcs;

  // ---------- setup ----------

  // create a buffer for storing temperature fields

  Vector A;
  int source_x;
  if (myrank == 0)
  {
    A = createVector(N);
    // set up initial conditions in A
    for (int i = 0; i < N; i++) {
      A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    source_x = N / 4;
    A[source_x] = 273 + 60;

    printf("Initial:\t");
    printTemperature(A, N);
    printf("\n");
  }
  MPI_Bcast(&A, N, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&source_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
  //value_t?
  //MPI_Scatter(A, elemsPerRank, MPI_INT, &A_, elemsPerRank, MPI_INT, 0, MPI_COMM_WORLD);

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(N);

  // send/receive requests for left/right rank
  MPI_Request LSrequest;
  MPI_Request LRrequest;
  MPI_Request RSrequest;
  MPI_Request RRrequest;

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (long long j = 0; j < elemsPerRank; j++)
    {
      int i = j + myrank * elemsPerRank;

      // center stays constant (the heat is still on)
      if (i == source_x)
      {
        B[i] = A[i];
        continue;
      }

      // get temperature at current position
      value_t tc = A[i];

      // get temperatures of adjacent cells
      value_t tl = (i != 0) ? A[i - 1] : tc;
      value_t tr = (i != N - 1) ? A[i + 1] : tc;

      // compute new temperature at current position
      B[i] = tc + 0.2 * (tl + tr + (-2 * tc));

      // send "data corners" of current rank to the prev and next rank respectively
      if (j == 0 && myrank != 0){
        MPI_Isend(&B[i], 1, MPI_INT, myrank - 1, 41, MPI_COMM_WORLD, &LSrequest);
        MPI_Irecv(&B[i-1], 1, MPI_INT, myrank - 1, 42, MPI_COMM_WORLD, &LRrequest);
      }
      else if (j == elemsPerRank-1 && myrank != numProcs-1){
        MPI_Isend(&B[i], 1, MPI_INT, myrank + 1, 42, MPI_COMM_WORLD, &RSrequest);
        MPI_Irecv(&B[i+1], 1, MPI_INT, myrank + 1, 41, MPI_COMM_WORLD, &RRrequest);
      }
    }

    // swap matrices (just pointers, not content)
    MPI_Barrier(MPI_COMM_WORLD);

    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (rank == 0 && !(t % 1000)) {
      printf("Step t=%d:\t", t);
      printTemperature(A, N);
      printf("\n");
    }
  }

  releaseVector(B);

  // ---------- check ----------

  printf("Final:\t\t");
  printTemperature(A, N);
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

  MPI_Finalize();

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int N) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;

  // step size in each dimension
  int sW = N / W;

  // room
  // left wall
  printf("X");
  // actual room
  for (int i = 0; i < W; i++) {
    // get max temperature in this tile
    value_t max_t = 0;
    for (int x = sW * i; x < sW * i + sW; x++) {
      max_t = (max_t < m[x]) ? m[x] : max_t;
    }
    value_t temp = max_t;

    // pick the 'color'
    int c = ((temp - min) / (max - min)) * numColors;
    c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

    // print the average temperature
    printf("%c", colors[c]);
  }
  // right wall
  printf("X");
}
