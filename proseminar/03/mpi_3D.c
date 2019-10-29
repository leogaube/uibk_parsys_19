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
  int dims[1] = {numProcs};
  int periods[1] = {0};
  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, &slices_2D);
  MPI_Comm_rank(slices_2D, &rank);
  printf("rank after init %d\n",rank);

  // get the adjacent slices
  int top_rank = rank;
  int bottom_rank = rank;
  MPI_Cart_shift(slices_2D, 0, 1, &bottom_rank, &top_rank);
  printf("rank after shift %d\n",rank);
  if(top_rank == MPI_PROC_NULL) { top_rank = rank; }
  if(bottom_rank == MPI_PROC_NULL) { bottom_rank = rank; }

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

  MPI_Isend(A, Nx*Ny, MPI_FLOAT, top_rank, 0, slices_2D, &topSRequest);
  MPI_Isend(&(A[IDX_3D(0,0,Mz-1,Nx,Ny)]), Nx*Ny, MPI_FLOAT, bottom_rank, 0, slices_2D, &bottomSRequest);

  MPI_Irecv(bottom_layer, Nx*Ny, MPI_FLOAT, bottom_rank, 0, slices_2D, &bottomRRequest);
  MPI_Irecv(top_layer, Nx*Ny, MPI_FLOAT, top_rank, 0, slices_2D, &topRRequest);

  // for each time step ..
  for (int t = 0; t < T; t++) {
	// send the uppermost and lowest layer to upper and lower slices
	// .. we propagate the temperature
    for (int z = 0; z < Mz; z++)
    {
	  //sync data
	  if (z == 0){
	    MPI_Wait(&bottomSRequest, MPI_STATUS_IGNORE);
	    MPI_Wait(&topRRequest, MPI_STATUS_IGNORE);
          } else if (z == Mz - 1){
	    MPI_Wait(&topSRequest, MPI_STATUS_IGNORE);
	    MPI_Wait(&bottomRRequest, MPI_STATUS_IGNORE);
  	  }
	  printf("wait done %d\n",rank);
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

      // send/receive "data corners" to/from the prev/next rank
      if (z == 0){
    	  MPI_Isend(A, Nx*Ny, MPI_FLOAT, top_rank, 0, slices_2D, &topSRequest);
    	  MPI_Irecv(bottom_layer, Nx*Ny, MPI_FLOAT, bottom_rank, 0, slices_2D, &bottomRRequest);
      }
      else if (z == Mz - 1){
    	  MPI_Isend(&(A[IDX_3D(0,0,Mz-1,Nx,Ny)]), Nx*Ny, MPI_FLOAT, bottom_rank, 0, slices_2D, &bottomSRequest);
    	  MPI_Irecv(top_layer, Nx*Ny, MPI_FLOAT, top_rank, 0, slices_2D, &topRRequest);
      }
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;
    printf("t/T: %d/%d %d\n",t,T,rank);
  }
  printf("DONE %d\n", rank);

  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);
  printf("ok? %d\n",rank);
  Vector AA = NULL;
  if(rank == top_rank){
	  printf("top is %d\n",rank);
	  AA = createVector(Nx*Ny*Nz);
  }
  printf("hey %d\n",rank);
  MPI_Gather(A, Nx*Ny*Mz, MPI_FLOAT, AA, Nx*Ny*Mz, MPI_FLOAT, 0, slices_2D);
  printf("final %d\n",rank);
  releaseVector(A);

  if (rank == top_rank){
#ifdef VERBOSE
		printf("Final:\t\t");
		printTemperature(AA, Nx, Ny, Nz);
		printf("\n");
#endif

	  // ---------- check ----------
	  double residual = is_verified_3D(AA, Nx, Ny, Nz, source_x, source_y, source_z, T);
	  printf("The maximal deviation from the 1D theory is %fK.", residual);

	  // ---------- cleanup ----------
	  releaseVector(AA);

    double end = MPI_Wtime();
    printf("The process took %g seconds to finish. \n", end - start);
  }

  MPI_Finalize();

  // done
  return EXIT_SUCCESS;
}
