#include <mpi.h>
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
  int T = MAX(Nx, Ny)*500;

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
  int dims[1] = {numProcs};
  int periods[1] = {0};
  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &stripes_1D);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef VERBOSE
  if (rank == 0)
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d for T=%d timesteps\n", Nx, Ny, T);
#endif

  // get the adjacent slices
  int top_rank = rank;
  int bottom_rank = rank;
  MPI_Cart_shift(stripes_1D, 0, 1, &top_rank, &bottom_rank);
  if(top_rank == MPI_PROC_NULL) { top_rank = rank; }
  if(bottom_rank == MPI_PROC_NULL) { bottom_rank = rank; }
  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx*My);
  Vector AA = NULL;
  if(rank == 0){
    AA = createVector(Nx*Ny);
  }

  // set up initial conditions in A
  for (long i = 0; i < Nx*My; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  if(rank*My <= source_y && (rank+1)*My > source_y){
    A[IDX_2D(source_x,source_y - rank*My,Nx)] = 273 + 60;
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx*My);
  Vector upper_layer = createVector(Nx);
  Vector lower_layer = createVector(Nx);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // send the uppermost and lowest layer to upper and lower slices
    MPI_Send(A, Nx, MPI_FLOAT, top_rank, 0, stripes_1D);
    MPI_Send(&(A[IDX_2D(0,My-1,Nx)]), Nx, MPI_FLOAT, bottom_rank, 0, stripes_1D);

    // .. we propagate the temperature
    for (int y = 0; y < My; y++)
    {
      if(y==0){
    	  MPI_Recv(upper_layer, Nx, MPI_FLOAT, top_rank, 0, stripes_1D, MPI_STATUS_IGNORE);
      }
      if (y==My-1){
    	  MPI_Recv(lower_layer, Nx, MPI_FLOAT, bottom_rank, 0, stripes_1D, MPI_STATUS_IGNORE);
      }
        for (int x = 0; x < Nx; x++)
        {
          // get the current idx
          long i = IDX_2D(x, y, Nx);

          // center stays constant (the heat is still on)
          if (i+rank*Nx*My == IDX_2D(source_x, source_y, Nx)) {
            B[i] = A[i];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i];

          // get temperatures of adjacent cells
          value_t tl = (x != 0) ? A[IDX_2D(x - 1, y, Nx)] : tc;
          value_t tr = (x != Nx - 1) ? A[IDX_2D(x + 1, y, Nx)] : tc;
	        value_t tu, td;
          if(rank == top_rank){
            tu = (y != 0) ? A[IDX_2D(x, y - 1, Nx)] : tc;
          } else {
            tu = (y != 0) ? A[IDX_2D(x, y - 1, Nx)] : upper_layer[x];
          }
          if(rank == bottom_rank){
            td = (y != My - 1) ? A[IDX_2D(x, y + 1, Nx)] : tc;
          } else {
            td = (y != My - 1) ? A[IDX_2D(x, y + 1, Nx)] : lower_layer[x];
          }

          // compute new temperature at current position
          B[i] = tc + 0.4/4 * (tl + tr + tu + td + (-4 * tc));
        }
      }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

#ifdef VERBOSE
    // show intermediate step
    if (!(t % 1000))
    {
      MPI_Gather(A, Nx*My, MPI_FLOAT, AA, Nx*My, MPI_FLOAT, 0, stripes_1D);
      if (rank == 0)
      {
        printf("Step t=%d:\n", t);
        printTemperature(AA, Nx, Ny, 1);
        printf("\n");
      }
    }
#endif
  }

  releaseVector(B);
  releaseVector(lower_layer);
  releaseVector(upper_layer);
  MPI_Gather(A, Nx*My, MPI_FLOAT, AA, Nx*My, MPI_FLOAT, 0, stripes_1D);
  releaseVector(A);

  if (rank == 0){
#ifdef VERBOSE
    printf("Final:\n");
    printTemperature(AA, Nx, Ny, 1);
    printf("\n");
#endif

    // ---------- check ----------
    double residual = is_verified_2D(AA, Nx, Ny, source_x, source_y, T);
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

