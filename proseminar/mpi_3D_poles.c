#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

void calculate_leftAreaToSend(value_t* src, value_t* dest, int Nx, int My, int Mz);

void calculate_rightAreaToSend(value_t* src, value_t* dest, int Nx, int My, int Mz);

void gatherCorrection(float* input, float* output, int numProcs, int Nx, int Ny, int Nz, int My, int Mz);

// -- simulation code ---

int main(int argc, char **argv)
{
  double start = MPI_Wtime();

  // 'parsing' optional input parameter = problem size
  int Nx = 10;
  int Ny = 10;
  int Nz = 10;
  if (argc == 2) {
    Nx = Ny = Nz = atoi(argv[1]);
  }
  int T = MAX(MAX(Nx, Ny), Nz) * 500;

  // MPI setup
  int rank, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (Nz % numProcs != 0)
  {
    printf("This Problem cannot be split up evenly among MPI ranks! (Nz mod numProcs != 0)");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  int Mz = Nz / sqrt(numProcs);
  int My = Ny / sqrt(numProcs);

  MPI_Comm poles_3D;
  int dims[2] = {0,0};
  MPI_Dims_create(numProcs,2,dims);
  int periods[2] = {0,0};
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &poles_3D);
  MPI_Comm_rank(poles_3D, &rank);

#ifdef VERBOSE
  if (rank == 0)
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
#endif

  // get the adjacent slices
  int top_rank = rank;
  int bottom_rank = rank;
  int left_rank = rank;
  int right_rank = rank;
  MPI_Cart_shift(poles_3D, 0, 1, &top_rank, &bottom_rank);
  MPI_Cart_shift(poles_3D, 1, 1, &left_rank, &right_rank);
  if (top_rank == MPI_PROC_NULL) {
    top_rank = rank;
  }   
  if (bottom_rank == MPI_PROC_NULL) {
    bottom_rank = rank;
  }
  if (left_rank == MPI_PROC_NULL) {
    left_rank = rank;
  }
  if (right_rank == MPI_PROC_NULL) {
    right_rank = rank;
  }

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx * My * Mz);

  // set up initial conditions in A
  for (long i = 0; i < Nx * My * Mz; i++)
  {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;

  if (rank * Mz <= source_z && (rank + 1) * Mz > source_z && rank * My <= source_y && (rank + 1) * My > source_y) {
    A[IDX_3D(source_x, source_y - (rank * My), source_z - (rank * Mz), Nx, My)] = 273 + 60;
  }


  // TODO: Test this stuff too
  Vector AA = NULL;
  Vector output = NULL;
  if (rank == 0)
  {
    AA = createVector(Nx * Ny * Nz);
    output = createVector(Nx * Ny * Nz);
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx * Ny * Mz);
  Vector top_layer = createVector(Nx * Ny);
  Vector bottom_layer = createVector(Nx * Ny);
  Vector left_layer = createVector(Nx * Ny);
  Vector right_layer = createVector(Nx * Ny);

  // for each time step ..
  for (int t = 0; t < T; t++)
  {
    MPI_Send(A, Nx*My, MPI_FLOAT, top_rank, 0, poles_3D);
    MPI_Send(&(A[IDX_3D(0,0,Mz-1,Nx,My)]), Nx*My, MPI_FLOAT, bottom_rank, 0, poles_3D);
		value_t *leftAreaToSend = createVector(Nx*Mz);
		calculate_leftAreaToSend(A,leftAreaToSend,Nx,My,Mz);
    MPI_Send(leftAreaToSend, Nx*Mz, MPI_FLOAT, left_rank, 0, poles_3D);
		value_t *rightAreaToSend = createVector(Nx*Mz);
		calculate_rightAreaToSend(A,rightAreaToSend,Nx,My,Mz);
    MPI_Send(rightAreaToSend, Nx*Mz, MPI_FLOAT, right_rank, 0, poles_3D);

		// If we use non blocking communication change this!! Releasing should then be done after wait
		releaseVector(leftAreaToSend);
		releaseVector(rightAreaToSend);

		MPI_Recv(top_layer, Nx*My, MPI_FLOAT, top_rank, 0, poles_3D, MPI_STATUS_IGNORE);
	  MPI_Recv(bottom_layer, Nx*My, MPI_FLOAT, bottom_rank, 0, poles_3D, MPI_STATUS_IGNORE);
	  MPI_Recv(left_layer, Nx*Mz, MPI_FLOAT, left_rank, 0, poles_3D, MPI_STATUS_IGNORE);
	  MPI_Recv(right_layer, Nx*Mz, MPI_FLOAT, right_rank, 0, poles_3D, MPI_STATUS_IGNORE);

    // send the uppermost and lowest layer to upper and lower slices
    // .. we propagate the temperature
    for (int z = 0; z < Mz; z++) {
      for (int y = 0; y < My; y++) {
        for (int x = 0; x < Nx; x++) {
          // get the current idx
          long i = IDX_3D(x, y, z, Nx, My);

          // center stays constant (the heat is still on)
          if (i == IDX_3D(source_x, source_y-(rank*My), source_z-(rank*Mz), Nx, My))
          {
            B[i] = A[i];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i];

          // get temperatures of adjacent cells
          value_t tl = (x != 0) ? A[IDX_3D(x - 1, y, z, Nx, Ny)] : tc;
          value_t tr = (x != Nx - 1) ? A[IDX_3D(x + 1, y, z, Nx, Ny)] : tc;
          value_t tu = (y != 0) ? A[IDX_3D(x, y - 1, z, Nx, Ny)] : left_layer[IDX_2D(x, y, Nx)];
          value_t td = (y != Ny - 1) ? A[IDX_3D(x, y + 1, z, Nx, Ny)] : right_layer[IDX_2D(x, y, Nx)];
          value_t tf = (z != 0) ? A[IDX_3D(x, y, z - 1, Nx, Ny)] : top_layer[IDX_2D(x, y, Nx)];
          value_t tb = (z != Mz - 1) ? A[IDX_3D(x, y, z + 1, Nx, Ny)] : bottom_layer[IDX_2D(x, y, Nx)];

          // compute new temperature at current position
          B[i] = tc + 0.16666 * (tl + tr + tu + td + tf + tb + (-6 * tc));
        }
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
      MPI_Gather(A, Nx * Ny * Mz, MPI_FLOAT, AA, Nx * Ny * Mz, MPI_FLOAT, 0, poles_3D);
      if (rank == 0)
      {
        printf("Step t=%d:\n", t);
				gatherCorrection(AA, output, numProcs, Nx, Ny, Nz, My, Mz);
        printTemperature(output, Nx, Ny, Nz);
        printf("\n");
      }
    }
#endif
  }

  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);
  releaseVector(left_layer);
  releaseVector(right_layer);

  MPI_Gather(A, Nx * Ny * Mz, MPI_FLOAT, AA, Nx * Ny * Mz, MPI_FLOAT, 0, poles_3D);
	gatherCorrection(AA, output, numProcs, Nx, Ny, Nz, My, Mz);
  releaseVector(A);

  if (rank == 0)
  {
#ifdef VERBOSE
    printf("Final:\n");
    printTemperature(output, Nx, Ny, Nz);
    printf("\n");
#endif

    // ---------- check ----------
    double residual = is_verified_3D(output, Nx, Ny, Nz, source_x, source_y, source_z, T);
    printf("The maximal deviation from the 1D theory is %fK.", residual);

    // ---------- cleanup ----------
    releaseVector(AA);
		releaseVector(output);

    double end = MPI_Wtime();
    printf("The process took %g seconds to finish. \n", end - start);
  }

  MPI_Finalize();

  // done
  return EXIT_SUCCESS;
}

void calculate_leftAreaToSend(value_t* src, value_t* dest, int Nx, int My, int Mz) {
  int k = 0;
	for(int i = 0; i < Mz; i ++) {
		for(int j = 0; j < Nx; j++) {
			dest[k] = src[IDX_3D(j,0,i,Nx,My)];
    	k++;
		}
  }
}

void calculate_rightAreaToSend(value_t* src, value_t* dest, int Nx, int My, int Mz) {
  int k = 0;
	for(int i = 0; i < Mz; i ++) {
		for(int j = 0; j < Nx; j++) {
			dest[k] = src[IDX_3D(j,My-1,i,Nx,My)];
    	k++;
		}
  }
}


void gatherCorrection(float* input, float* output, int numProcs, int Nx, int Ny, int Nz, int My, int Mz) {
	int l = 0;  
	for(int h = 0; h < numProcs; h++) {
		for(int i = 0; i < Mz; i++) {
		  for(int j = 0; j < My; j++) {
		    for(int k = 0; k < Nx; k++) {
					output[IDX_3D(k, j+(h*My), i+(h*My), Nx, Ny)] = input[l];
					l++;
		    }
		  }
		}
	}
}
