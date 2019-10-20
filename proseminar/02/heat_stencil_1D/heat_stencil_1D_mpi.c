#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

typedef double value_t;

#define RESOLUTION 120

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

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
  int T = 10000;

  // MPI setup
  int rank, numProcs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (N % numProcs != 0){
    printf("This Problem cannot be split up evenly among MPI ranks! (N mod numProcs != 0)");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  int M = N / numProcs;

  if (rank == 0)
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps using %d processes with subroom size M=%d\n", N, T, numProcs, M);

  // ---------- setup ----------

  // create a buffer for storing temperature fields

  Vector AA;
  Vector A = createVector(M);
  int source_x;
  if (rank == 0)
  {
    AA = createVector(N);
    // set up initial conditions in A
    for (int i = 0; i < N; i++) {
      AA[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    source_x = N/4;
    AA[source_x] = 273 + 60;

    printf("Initial:\t");
    printTemperature(AA, N);
    printf("\n");
  }
  MPI_Scatter(AA, M, MPI_DOUBLE, A, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&source_x, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(M);

  // send/receive requests for left/right rank
  /*
  MPI_Request LSrequest;
  MPI_Request LRrequest;
  MPI_Request RSrequest;
  MPI_Request RRrequest;
  */
  MPI_Request LRrequest;
  MPI_Request RRrequest;

  double leftCell;
  double rightCell;

  MPI_Bsend(&(A[0]), 1, MPI_DOUBLE, MAX(rank-1, 0), 0, MPI_COMM_WORLD);
  MPI_Bsend(&(A[M - 1]), 1, MPI_DOUBLE, MIN(rank+1, numProcs-1), 0, MPI_COMM_WORLD);

  MPI_Irecv(&leftCell, 1, MPI_DOUBLE, MAX(rank - 1, 0), 0, MPI_COMM_WORLD, &LRrequest);
  MPI_Irecv(&rightCell, 1, MPI_DOUBLE, MIN(rank + 1, numProcs - 1), 0, MPI_COMM_WORLD, &RRrequest);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (long long i = 0; i < M; i++)
    {
      // center stays constant (the heat is still on)
      if (i + (rank * M) == source_x)
      {
        B[i] = A[i];
      }
      else{
        //sync data
        if (i == 0)
          MPI_Wait(&LRrequest, MPI_STATUS_IGNORE);
        else if (i == M - 1)
          MPI_Wait(&RRrequest, MPI_STATUS_IGNORE);

        // get temperature at current position
        value_t tc = A[i];

        // get temperatures of adjacent cells
        value_t tl = (i != 0) ? A[i - 1] : leftCell;
        value_t tr = (i != M - 1) ? A[i + 1] : rightCell;

        // compute new temperature at current position
        B[i] = tc + 0.2 * (tl + tr + (-2 * tc));
      }

      // send/receive "data corners" to/from the prev/next rank
      if (i == 0){
        MPI_Bsend(&(B[i]), 1, MPI_DOUBLE, MAX(rank - 1, 0), 0, MPI_COMM_WORLD);
        MPI_Irecv(&leftCell, 1, MPI_DOUBLE, MAX(rank - 1, 0), 0, MPI_COMM_WORLD, &LRrequest);
      }
      else if (i == M-1){
        MPI_Bsend(&(B[i]), 1, MPI_DOUBLE, MIN(rank + 1, numProcs - 1), 0, MPI_COMM_WORLD);
        MPI_Irecv(&rightCell, 1, MPI_DOUBLE, MIN(rank + 1, numProcs - 1), 0, MPI_COMM_WORLD, &RRrequest);
      }
      //printf("timestep: %d, %lld, %d\n", t, i, rank);
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 1))
    {
      MPI_Gather(A, M, MPI_DOUBLE, AA, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (rank == 0)
      {
        printf("Step t=%d:\t", t);
        printTemperature(A, M);
        printf("\n");
      }
    }
  }

  releaseVector(B);

  // ---------- check ----------

  MPI_Gather(A, M, MPI_DOUBLE, AA, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0){
    printf("Final:\t\t");
    printTemperature(AA, N);
    printf("\n");
  }

  int success = 1;
  for (long long i = 0; i < M; i++) {
    value_t temp = A[i];
    if (273 <= temp && temp <= 273 + 60)
      continue;
    success = 0;
    break;
  }

  printf("Verification for rank %d: %s\n", rank, (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  if (rank == 0)
  {
    releaseVector(AA);
  }
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
