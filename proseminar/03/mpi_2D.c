#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>

#include "heat_stencil.h"

#define RESOLUTION 120

void printTemperature(Vector m, int nx, int ny);

void calculate_leftCellsToSend(double* src, double* dest, int M);

void calculate_rightCellsToSend(double* src, double* dest, int M);


// -- simulation code ---

int main(int argc, cha#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

// -- simulation code ---

int main(int argc, char **argv) {
  double start = MPI_Wtime();

  // 'parsing' optional input parameter = problem size
  int Nx = 10;
  int Ny = 10;
  if (argc == 2) {
    Nx = Ny = atoi(argv[1]);
  }
  int T = MAX(Nx, Ny)*100;
#ifdef VERBOSE
  printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
#endif

  // MPI setup
  int rank, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (Ny % numProcs != 0){
    printf("This Problem cannot be split up evenly among MPI ranks! (Nz mod numProcs != 0)");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  int My = Ny / numProcs;

  MPI_Comm stripes_1D;
  int dims[1] = {1};
  int periods[1] = {0};
  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &stripes_1D);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // get the adjacent slices
  int top_rank = rank;
  int bottom_rank = rank;
  MPI_Cart_shift(slices_2D, 0, 1, &bottom_rank, &top_rank);
  if(top_rank == MPI_PROC_NULL) { top_rank = rank; }
  if(bottom_rank == MPI_PROC_NULL) { bottom_rank = rank; }

  MPI_Request topSRequest;
  MPI_Request topRRequest;
  MPI_Request bottomSRequest;
  MPI_Request bottomRRequest;
  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx*My);

  // set up initial conditions in A
  for (long i = 0; i < Nx*My; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  if(rank*My<=source_y && (rank+1)*My>source_y){
	  A[IDX_2D(source_x,source_y,Ny)] = 273 + 60;
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx*My);
  Vector upper_layer = createVector(Nx);
  Vector lower_layer = createVector(Nx);

  // for each time step ..
  for (int t = 0; t < T; t++) {

	// send the uppermost and lowest layer to upper and lower slices
	MPI_Isend(A, Nx, MPI_FLOAT, top_rank, 0, stripes_1D, &topSRequest);
	MPI_Isend(&(A[IDX_2D(0,My-1,Nx)]), Nx, MPI_FLOAT, bottom_rank, 0, stripes_1D, &bottomSRequest);

	// .. we propagate the temperature
    for (int y = 0; y < My; y++)
    {
      if(y==0){
    	  MPI_Irecv(upper_layer, Nx, MPI_FLOAT, top_rank, 0, stripes_1D, &topRRequest);
      } else if (y==My-1){
    	  MPI_Irecv(lower_layer, Nx, MPI_FLOAT, bottom_rank, 0, stripes_1D, &bottomRRequest);
      }
        for (int x = 0; x < Nx; x++)
        {
          // get the current idx
          long i = IDX_2D(x, y, Nx);

          // center stays constant (the heat is still on)
          if (i == IDX_2D(source_x, source_y, Nx))
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
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

  }

  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);
  Vector AA = NULL;
  if(rank == 0){
	  AA = createVector(Nx*Ny*Nz);
  }
  MPI_Gather(A, Nx*Ny*Mz, MPI_FLOAT, AA, Nx*Ny*Mz, MPI_FLOAT, 0, slices_2D);
  releaseVector(A);

  if (rank == 0){
#ifdef VERBOSE
		printf("Final:\t\t");
		printTemperature(AA, Nx, Ny, Nz);
		printf("\n");
#endif

	  // ---------- check ----------
	  double residual = is_verified_3D(AA, Nx, Ny, Nz, source_x, source_y, source_z, T);
	  printf("The maximal deviation from the 1D theory is %fK.", residual);
  }

  // ---------- cleanup ----------
  releaseVector(AA);

  if (rank == 0)
  {
    double end = MPI_Wtime();
    printf("The process took %g seconds to finish. \n", end - start);
  }

  MPI_Finalize();

  // done
  return EXIT_SUCCESS;
}
r **argv) {
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
    double #include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

// -- simulation code ---

int main(int argc, char **argv) {
  double start = MPI_Wtime();

  // 'parsing' optional input parameter = problem size
  int Nx = 10;
  int Ny = 10;
  int Nz = 10;
  if (argc == 2) {
    Nx = Ny = Nz = atoi(argv[1]);
  } else if (argc == 4){
	  Nx = atoi(argv[1]);
	  Ny = atoi(argv[2]);
	  Nz = atoi(argv[3]);
  }
  int T = MAX(MAX(Nx, Ny), Nz)*100;
#ifdef VERBOSE
  printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
#endif

  // MPI setup
  int rank, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (Nz % numProcs != 0){
    printf("This Problem cannot be split up evenly among MPI ranks! (Nz mod numProcs != 0)");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  int Mz = Nz / numProcs;

  MPI_Comm slices_2D;
  int dims[1] = {1};
  int periods[1] = {0};
  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &slices_2D);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // get the adjacent slices
  int top_rank = rank;
  int bottom_rank = rank;
  if(rank>0){
	  MPI_Cart_shift(slices_2D, 0, 1, &rank, &top_rank);
  }
  if(rank<numProcs-1){
	  MPI_Cart_shift(slices_2D, 0, -1, &rank, &bottom_rank);
  }

  MPI_Request topSRequest;
  MPI_Request topRRequest;
  MPI_Request bottomSRequest;
  MPI_Request bottomRRequest;
  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx*Ny*Mz);

  // set up initial conditions in A
  for (long i = 0; i < Nx*Ny*Mz; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;
  if(rank*Mz<=source_z && (rank+1)*Mz>source_z){
	  A[IDX_3D(source_x,source_y, source_z,Nx,Ny)] = 273 + 60;
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx*Ny*Mz);
  Vector top_layer = createVector(Nx*Ny);
  Vector bottom_layer = createVector(Nx*Ny);

  // for each time step ..
  for (int t = 0; t < T; t++) {

	// send the uppermost and lowest layer to upper and lower slices
	MPI_Isend(A, Nx*Ny, MPI_FLOAT, top_rank, 0, slices_2D, &topSRequest);
	MPI_Isend(&(A[IDX_3D(0,0,Mz-1,Nx,Ny)]), Nx*Ny, MPI_FLOAT, bottom_rank, 0, slices_2D, &bottomSRequest);

	// .. we propagate the temperature
    for (int z = 0; z < Mz; z++)
    {
      if(z==0){
    	  MPI_Irecv(top_layer, Nx*Ny, MPI_FLOAT, top_rank, 0, slices_2D, &topRRequest);
      } else if (z==Mz-1){
    	  MPI_Irecv(bottom_layer, Nx*Ny, MPI_FLOAT, bottom_rank, 0, slices_2D, &bottomRRequest);
      }
      for (int y = 0; y < Ny; y++)
      {
        for (int x = 0; x < Nx; x++)
        {
          // get the current idx
          long i = IDX_3D(x, y, z, Nx, Ny);

          // center stays constant (the heat is still on)
          if (i == IDX_3D(source_x, source_y, source_z, Nx, Ny))
          {
            B[i] = A[i];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i];

          // get temperatures of adjacent cells
          value_t tl = (x != 0) ? A[IDX_3D(x - 1, y, z, Nx, Ny)] : tc;
          value_t tr = (x != Nx - 1) ? A[IDX_3D(x + 1, y, z, Nx, Ny)] : tc;
          value_t tu = (y != 0) ? A[IDX_3D(x, y - 1, z, Nx, Ny)] : tc;
          value_t td = (y != Ny - 1) ? A[IDX_3D(x, y + 1, z, Nx, Ny)] : tc;
          value_t tf = (z != 0) ? A[IDX_3D(x, y, z - 1, Nx, Ny)] : top_layer[IDX_2D(x, y, Nx)];
          value_t tb = (z != Mz - 1) ? A[IDX_3D(x, y, z + 1, Nx, Ny)] : bottom_layer[IDX_2D(x, y, Nx)];

          // compute new temperature at current position
          B[i] = tc + 0.4/6 * (tl + tr + tu + td + tf + tb + (-6 * tc));
        }
      }
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

  }

  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);
  Vector AA = NULL;
  if(rank == 0){
	  AA = createVector(Nx*Ny*Nz);
  }
  MPI_Gather(A, Nx*Ny*Mz, MPI_FLOAT, AA, Nx*Ny*Mz, MPI_FLOAT, 0, slices_2D);
  releaseVector(A);

  if (rank == 0){
#ifdef VERBOSE
		printf("Final:\t\t");
		printTemperature(AA, Nx, Ny, Nz);
		printf("\n");
#endif

	  // ---------- check ----------
	  double residual = is_verified_3D(AA, Nx, Ny, Nz, source_x, source_y, source_z, T);
	  printf("The maximal deviation from the 1D theory is %fK.", residual);
  }

  // ---------- cleanup ----------
  releaseVector(AA);

  if (rank == 0)
  {
    double end = MPI_Wtime();
    printf("The process took %g seconds to finish. \n", end - start);
  }

  MPI_Finalize();

  // done
  return EXIT_SUCCESS;
}
leftCellsToSend[M];
    calculate_leftCellsToSend(A, leftCellsToSend, M);
    MPI_Isend(&(leftCellsToSend[0]), M, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &LSrequest);
  }
  if (rank % numberOfRanksPerRow != numberOfRanksPerRow - 1) {
    double rightCellsToSend[M];
    calculate_rightCellsToSend(A, rightCellsToSend, M);
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
          value_t tr = rank % numberOfRanksPerRow != numberOfRanksPerRow - 1 ? (x != M - 1) ? A[IDX_2D(x+1,y,N)] : rightCells[y] : A[i];
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
      calculate_leftCellsToSend(A, leftCellsToSend, M);
      MPI_Isend(&(leftCellsToSend[0]), M, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &LSrequest);
    }
    if (rank % numberOfRanksPerRow != numberOfRanksPerRow - 1) {
      double rightCellsToSend[M];
      calculate_rightCellsToSend(A, rightCellsToSend, M);
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
        printTemperature(AA, N, N);
        printf("\n");
      }
    }
  }

  releaseVector(B);

  // ---------- check ----------

  MPI_Gather(A, M, MPI_DOUBLE, AA, M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (rank == 0){
    printf("Final:\t\t");
    printTemperature(AA, N, N);
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
