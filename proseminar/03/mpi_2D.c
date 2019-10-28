#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>

#include <time.h>

typedef double value_t;

#define RESOLUTION 120

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int N);

void calculate_leftCellsToSend(double* src, double* dest, int M);

void calculate_rightCellsToSend(double* src, double* dest, int M);
// -- simulation code ---

int main(int argc, char **argv) {
  double start = MPI_Wtime();

  // 'parsing' optional input parameter = problem size
  int N = 128;
  if (argc > 1)
  {
    N = atoi(argv[1]);
  }
  int T = N * 500;

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

  int M = sqrt(N*N / numProcs);

  if (rank == 0)
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps using %d processes with subroom size M=%d\n", N, T, numProcs, M);

  // ---------- setup ----------

  // create a buffer for storing temperature fields

  Vector AA = NULL;
  Vector A = createVector(M * M);
  int source_x;
  int source_y;
  if (rank == 0)
  {
    AA = createVector(N*N);
    // set up initial conditions in A
    for (int i = 0; i < N*N; i++) {
      AA[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source in one corner
    source_x = N*N / 4;
    source_y = source_x;
    AA[IDX_2D(source_x,source_y,N*N)] = 273 + 60;

#ifdef VERBOSE
  	printf("time %i:\n",t);
	  printTemperature(B, N, N);
#endif
  }
  MPI_Scatter(AA, M*M, MPI_DOUBLE, A, M*M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&source_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&source_y, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(M*M);

  // send/receive requests for left/right/top/bottom rank
  MPI_Request LRrequest;
  MPI_Request LSrequest;
  MPI_Request RRrequest;
  MPI_Request RSrequest;
  MPI_Request TRrequest;
  MPI_Request TSrequest;
  MPI_Request DRrequest;
  MPI_Request DSrequest;

  double leftCells[M];
  double rightCells[M];
  double topCells[M];
  double bottomCells[M];

  int numberOfRanksPerRow = sqrt(numProcs);

  if (rank % numberOfRanksPerRow != 0) {
    double leftCellsToSend[M];
    calculate_leftCellsToSend(A, leftCellsToSend, M)
    MPI_Isend(&(leftCellsToSend[0]), M, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &LSrequest);
  }
  if (rank % numberOfRanksPerRow != numberOfRanksPerRow - 1) {
    double rightCellsToSend[M];
    calculate_rightCellsToSend(A, rightCellsToSend, M)
    MPI_Isend(&(rightCellsToSend[0]), M, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &RSrequest);
  }
  if (rank - numberOfRanksPerRow > 0) {
    MPI_Isend(&(A[0]), M, MPI_DOUBLE, rank - numberOfRanksPerRow, 0, MPI_COMM_WORLD, &TSrequest);
  }
  if (rank + numberOfRanksPerRow < numberOfRanksPerRow * numberOfRanksPerRow) {
    MPI_Isend(&(A[M*M-M-1]), M, MPI_DOUBLE, rank + numberOfRanksPerRow, 0, MPI_COMM_WORLD, &DSrequest);
  }
  
  if (rank % numberOfRanksPerRow != 0) {
    MPI_Irecv(&leftCells, M, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, &LRrequest);
  }
  if (rank % numberOfRanksPerRow != numberOfRanksPerRow - 1) {
    MPI_Irecv(&rightCells, M, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &RRrequest);
  }
  if (rank - numberOfRanksPerRow > 0) {
    MPI_Irecv(&topCells, M, MPI_DOUBLE, rank - numberOfRanksPerRow, 0, MPI_COMM_WORLD, &TRrequest);
  }
  if (rank + numberOfRanksPerRow < numberOfRanksPerRow * numberOfRanksPerRow) {
    MPI_Irecv(&bottomCells, M, MPI_DOUBLE, rank + numberOfRanksPerRow, 0, MPI_COMM_WORLD, &DRrequest);
  }


  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    //sync data
    MPI_Wait(&LSrequest, MPI_STATUS_IGNORE);
    MPI_Wait(&LRrequest, MPI_STATUS_IGNORE);
  
    MPI_Wait(&RSrequest, MPI_STATUS_IGNORE);
    MPI_Wait(&RRrequest, MPI_STATUS_IGNORE);

    MPI_Wait(&TSrequest, MPI_STATUS_IGNORE);
    MPI_Wait(&TRrequest, MPI_STATUS_IGNORE);

    MPI_Wait(&DSrequest, MPI_STATUS_IGNORE);
    MPI_Wait(&DRrequest, MPI_STATUS_IGNORE);

    for (int y = 0; y < M; y++) {
      for (int x = 0; x < M; x++) {
        int i = IDX_2D(x,y,N);
        // center stays constant (the heat is still on)
        if (x + (rank * M) == source_x && y + (rank * M) == source_y) {
          B[i] = A[i];
        }
        else {
          // get temperature at current position
          value_t tc = A[i];

          // get temperatures of adjacent cells
          value_t tl = rank % numberOfRanksPerRow != 0 ? (x != 0) ? A[IDX_2D(x-1,y,N)] : leftCells[y] : A[i];
          value_t tr = rank % numberOfRanksPerRow != numberOfRanksPerRow - 1) ? (x != M - 1) ? A[IDX_2D(x+1,y,N)] : rightCells[y] : A[i];
          value_t tu = rank - numberOfRanksPerRow > 0 ? (y != 0) ? A[IDX_2D(x,y-1,N)] : topCells[x] : A[i];
          value_t tb = rank + numberOfRanksPerRow < numberOfRanksPerRow * numberOfRanksPerRow ? (y != M - 1) ? A[IDX_2D(x,y+1,N)] : bottomCells[x] : A[i];

          // compute new temperature at current position
		    	B[i] = tc + 0.2 * (tl + tr + tu + tb + (-4 * tc));
        }
      }
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    if (rank % numberOfRanksPerRow != 0) {
      double leftCellsToSend[M];
      calculate_leftCellsToSend(A, leftCellsToSend, M)
      MPI_Isend(&(leftCellsToSend[0]), M, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &LSrequest);
    }
    if (rank % numberOfRanksPerRow != numberOfRanksPerRow - 1) {
      double rightCellsToSend[M];
      calculate_rightCellsToSend(A, rightCellsToSend, M)
      MPI_Isend(&(rightCellsToSend[0]), M, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &RSrequest);
    }
    if (rank - numberOfRanksPerRow > 0) {
      MPI_Isend(&(A[0]), M, MPI_DOUBLE, rank - numberOfRanksPerRow, 0, MPI_COMM_WORLD, &TSrequest);
    }
    if (rank + numberOfRanksPerRow < numberOfRanksPerRow * numberOfRanksPerRow) {
      MPI_Isend(&(A[M*M-M-1]), M, MPI_DOUBLE, rank + numberOfRanksPerRow, 0, MPI_COMM_WORLD, &DSrequest);
    }
    
    if (rank % numberOfRanksPerRow != 0) {
      MPI_Irecv(&leftCells, M, MPI_DOUBLE, rank -1, 0, MPI_COMM_WORLD, &LRrequest);
    }
    if (rank % numberOfRanksPerRow != numberOfRanksPerRow - 1) {
      MPI_Irecv(&rightCells, M, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &RRrequest);
    }
    if (rank - numberOfRanksPerRow > 0) {
      MPI_Irecv(&topCells, M, MPI_DOUBLE, rank - numberOfRanksPerRow, 0, MPI_COMM_WORLD, &TRrequest);
    }
    if (rank + numberOfRanksPerRow < numberOfRanksPerRow * numberOfRanksPerRow) {
      MPI_Irecv(&bottomCells, M, MPI_DOUBLE, rank + numberOfRanksPerRow, 0, MPI_COMM_WORLD, &DRrequest);
    }

    // show intermediate step
    if (!(t % 1000))
    {
      MPI_Gather(A, M, MPI_DOUBLE, AA, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (rank == 0)
      {
        printf("Step t=%d:\t", t);
        printTemperature(AA, N);
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

  MPI_Barrier(MPI_COMM_WORLD);

  int success = (is_verified_2D(A, N, N, source_x, source_y, T)==0);
  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  if (rank == 0)
  {
    releaseVector(AA);
  }
  releaseVector(A);

  if (rank == 0)
  {
    double end = MPI_Wtime();
    printf("The process took %g seconds to finish. \n", end - start);
  }

  MPI_Finalize();

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
		printf("%c", colors[c]);
	  }
	  // right wall
	  printf("X\n");
  }
}

void calculate_leftCellsToSend(double* src, double* dest, int M) {
  int k = 0;
  for(int i = 0; i < M*M; i += M) {
    dest[k] = src[i];
    k++;
  }
}

void calculate_rightCellsToSend(double* src, double* dest, int M) {
  int k = 0;
  for(int i = M-1; i < M*M; i +=M) {
    dest[k] = src[i];
    k++;
  }
}