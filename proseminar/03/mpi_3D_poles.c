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

  // get the adjacent slices
  int top_rank = rank;
  int bottom_rank = rank;
  int front_rank = rank;
  int back_rank = rank;
  MPI_Cart_shift(poles_3D, 0, 1, &top_rank, &bottom_rank);
  MPI_Cart_shift(poles_3D, 1, 1, &front_rank, &back_rank);
  if (top_rank == MPI_PROC_NULL) {
    top_rank = rank;
  }   
  if (bottom_rank == MPI_PROC_NULL) {
    bottom_rank = rank;
  }
  if (front_rank == MPI_PROC_NULL) {
    front_rank = rank;
  }
  if (back_rank == MPI_PROC_NULL) {
    back_rank = rank;
  }

	// 2D subroom slices for each dimension --> useful for sending 2D ghost cell slices to neighbouring ranks
	MPI_Datatype y_slice, z_slice;
  MPI_Type_vector(Mz, Nx, Nx * My, MPI_FLOAT, &y_slice);
  MPI_Type_vector(1, Nx * My, 1, MPI_FLOAT, &z_slice);
  MPI_Type_commit(&y_slice);
  MPI_Type_commit(&z_slice);

	// create subroom_datatypes for MPI_Gatherv
  //--> start of subroom determined by displacement_array
  //--> end of subroom determined by resized_subroom datatype
  int room_sizes[3] = {Nx, Ny, Nz};
  int subroom_sizes[3] = {Nx, My, Mz};
  int start_array[3] = {0, 0, 0};
  MPI_Datatype pole_subroom, resized_subroom;
  MPI_Type_create_subarray(3, room_sizes, subroom_sizes, start_array, MPI_ORDER_C, MPI_FLOAT, &pole_subroom);
  // 'pretend' that subroom is only 1 floats in size
  MPI_Type_create_resized(pole_subroom, 0, 1 * sizeof(float), &resized_subroom);
  MPI_Type_commit(&resized_subroom);

	int* recv_count_array = malloc(sizeof(int) * numProcs);
  int*	displacement_array = malloc(sizeof(int) * numProcs);
	if (rank == 0) {
		for (int r=0; r<numProcs; r++){
		  recv_count_array[r] = 1;
		  // get the global index of every rank at its local position 0
		  displacement_array[r] = local2global(r, 0, Nx, My, Mz, 1, sqrt(numProcs));
			//printf("local2global value: %d\n", local2global(r, 0, Nx, My, Mz, 1, sqrt(numProcs)));
		}
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

  Vector AA = NULL;
  if (rank == 0) {
    AA = createVector(Nx * Ny * Nz);
  }

	#ifdef VERBOSE
  MPI_Gatherv(A, Nx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, poles_3D);
  if (rank == 0) {
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);

    printf("Initial:\n");
    printTemperature(AA, Nx, Ny, Nz);
    printf("\n");
  }
#endif

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx * Ny * Mz);
  Vector top_layer = createVector(Nx * My);
  Vector bottom_layer = createVector(Nx * My);
  Vector front_layer = createVector(Nx * Mz);
  Vector back_layer = createVector(Nx * Mz);

  // for each time step ..
  for (int t = 0; t < T; t++)
  {
    MPI_Send(A, 1, z_slice, top_rank, 0, poles_3D);
    MPI_Send(&(A[IDX_3D(0,0,Mz-1,Nx,My)]), 1, z_slice, bottom_rank, 0, poles_3D);
    MPI_Send(A, 1, y_slice, front_rank, 0, poles_3D);
    MPI_Send(&(A[IDX_3D(0,My-1,0,Nx,My)]), 1, y_slice, back_rank, 0, poles_3D);

		MPI_Recv(top_layer, Nx*My, MPI_FLOAT, top_rank, 0, poles_3D, MPI_STATUS_IGNORE);
	  MPI_Recv(bottom_layer, Nx*My, MPI_FLOAT, bottom_rank, 0, poles_3D, MPI_STATUS_IGNORE);
	  MPI_Recv(front_layer, Nx*Mz, MPI_FLOAT, front_rank, 0, poles_3D, MPI_STATUS_IGNORE);
	  MPI_Recv(back_layer, Nx*Mz, MPI_FLOAT, back_rank, 0, poles_3D, MPI_STATUS_IGNORE);

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
          value_t tl = (x != 0) ? A[IDX_3D(x - 1, y, z, Nx, My)] : tc;
          value_t tr = (x != Nx - 1) ? A[IDX_3D(x + 1, y, z, Nx, My)] : tc;
        
          value_t tu = (z != 0) ? A[IDX_3D(x, y, z - 1, Nx, My)] : top_layer[IDX_2D(x, y, Nx)];
          value_t td = (z != Mz - 1) ? A[IDX_3D(x, y, z + 1, Nx, My)] : bottom_layer[IDX_2D(x, y, Nx)];
					value_t tf = (y != 0) ? A[IDX_3D(x, y - 1, z, Nx, My)] : front_layer[IDX_2D(x, z, Nx)];
          value_t tb = (y != My - 1) ? A[IDX_3D(x, y + 1, z, Nx, My)] : back_layer[IDX_2D(x, z, Nx)];

          // compute new temperature at current position
          B[i] = tc + 0.16666 * (tl + tr + tu + td + tf + tb + (-6 * tc));
        }
      }
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

		// show intermediate step
    if (!(t % 1000)) {
      MPI_Gatherv(A, Nx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, poles_3D);

      if (rank == 0) {
        printf("t: %d\n", t);
        printTemperature(AA, Nx, Ny, Nz);
        printf("\n");
      }
    }
  }

  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);
  releaseVector(front_layer);
  releaseVector(back_layer);

	//printf("rank :%d\n", rank);
  //printTemperature(A, Nx, My, Mz);
  //printf("\n");

  MPI_Gather(A, Nx * My * Mz, MPI_FLOAT, AA, Nx * Ny * Nz, MPI_FLOAT, 0, poles_3D);
	free(recv_count_array);
	free(displacement_array);
	
  releaseVector(A);

  if (rank == 0)
  {
/*
#ifdef VERBOSE
    printf("Final:\n");
    printTemperature(AA, Nx, Ny, Nz);
    printf("\n");
#endif


    // ---------- check ----------
    double residual = is_verified_3D(AA, Nx, Ny, Nz, source_x, source_y, source_z, T);
    printf("The maximal deviation from the 1D theory is %fK.", residual);
*/
    // ---------- cleanup ----------
    releaseVector(AA);

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
