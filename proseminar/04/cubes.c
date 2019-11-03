#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

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
  } else if (argc == 4) {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
  }
  int T = MAX(MAX(Nx, Ny), Nz) * TIMESTEPS_MUL;

  // MPI setup
  int rank, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  double cbrt_num_procs = cbrt((double)numProcs);
  int C = ceilf(cbrt_num_procs);
  if (C != (int)cbrt_num_procs) {
    printf("The number von ranks cannot be split up evenly into cubic subrooms! (numProcs != 8 || numProcs != 27 || numProcs != 64)");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  if (Nx % C != 0 || Ny % C != 0 || Nz % C != 0) {
    printf("This Problem cannot be split up evenly among MPI ranks in a cubic domain! (Nx mod C != 0 || Ny mod C != 0 || Nz mod C != 0)");
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  //TODO different values for Cx, Cy, Cz
  int Cx = C;
  int Cy = C;
  int Cz = C;

  int Mx = Nx / Cx;
  int My = Ny / Cy;
  int Mz = Nz / Cz;


  MPI_Request topRrequest, topSrequest, bottomRrequest, bottomSrequest;
  MPI_Request frontRrequest, frontSrequest, backRrequest, backSrequest;
  MPI_Request leftRrequest, leftSrequest, rightRrequest, rightSrequest;

  MPI_Comm cubes;

  int dims[3] = {Cx, Cy, Cz};
  int periods[3] = {0, 0, 0};

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cubes);
  MPI_Comm_rank(cubes, &rank);

#ifdef VERBOSE
  if (rank == 0){
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
    printf("C: %d, Mx: %d\n\n", C, Mx);
  }
#endif

  // get the adjacent slices
  int left_rank = rank;
  int right_rank = rank;
  int front_rank = rank;
  int back_rank = rank;
  int top_rank = rank;
  int bottom_rank = rank;

  MPI_Cart_shift(cubes, 0, 1, &top_rank, &bottom_rank);
  MPI_Cart_shift(cubes, 1, 1, &front_rank, &back_rank);
  MPI_Cart_shift(cubes, 2, 1, &left_rank, &right_rank);

  // possibly send data to yourself
  if (left_rank == MPI_PROC_NULL)
    left_rank = rank;
  if (right_rank == MPI_PROC_NULL)
    right_rank = rank;
  if (front_rank == MPI_PROC_NULL)
    front_rank = rank;
  if (back_rank == MPI_PROC_NULL)
    back_rank = rank;
  if (top_rank == MPI_PROC_NULL)
    top_rank = rank;
  if (bottom_rank == MPI_PROC_NULL)
    bottom_rank = rank;

#ifdef DEBUG
  printf("me: %d, l: %d, r: %d, f: %d, b: %d, u: %d, d: %d\n",
         rank, left_rank, right_rank, front_rank, back_rank, top_rank, bottom_rank);
#endif

  MPI_Datatype x_slice, y_slice, z_slice;
  MPI_Type_vector(My * Mz, 1, Mx, MPI_FLOAT, &x_slice);
  MPI_Type_vector(Mz, Mx, Mx * My, MPI_FLOAT, &y_slice);
  MPI_Type_vector(1, Mx * My, 1, MPI_FLOAT, &z_slice);
  MPI_Type_commit(&x_slice);
  MPI_Type_commit(&y_slice);
  MPI_Type_commit(&z_slice);


  //--> start of subroom determined by displacement_array
  //--> end of subroom determined by resized_subroom datatype
  int room_sizes[3] = {Nx, Ny, Nz};
  int subroom_sizes[3] = {Mx, My, Mz};
  int start_array[3] = {0, 0, 0};
  MPI_Datatype cubic_subroom, resized_subroom;
  MPI_Type_create_subarray(3, room_sizes, subroom_sizes, start_array, MPI_ORDER_C, MPI_FLOAT, &cubic_subroom);
  // 'pretend' that subroom is only Mx floats in size
  MPI_Type_create_resized(cubic_subroom, 0, 1 * sizeof(float), &resized_subroom);
  MPI_Type_commit(&resized_subroom);

  int* recv_count_array = malloc(sizeof(int) * numProcs);
  int* displacement_array = malloc(sizeof(int) * numProcs);
  for (int r=0; r<numProcs; r++){
    recv_count_array[r] = 1;
    // get the global index of every rank at its local position 0
    displacement_array[r] = local2global(r, 0, Mx, My, Mz, Cx, Cy);
  }

  // ---------- setup ----------

  // create a buffer for storing temperature fields for sub-rooms
  Vector A = createVector(Mx * My * Mz);

  // set up initial conditions in A
  for (int i = 0; i < Mx * My * Mz; i++) {
    A[i] = 273;  // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;
  int global_source = IDX_3D(source_x, source_y, source_z, Nx, Ny);
  int local_source = IDX_3D(source_x % Mx, source_y % My, source_z % Mz, Mx, My);

  if (local2global(rank, local_source, Mx, My, Mz, Cx, Cy) == global_source) {
    A[local_source] = 273 + 60;
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Mx * My * Mz);

  // create ghost exchange layers
  Vector left_layer = createVector(My * Mz);
  Vector right_layer = createVector(My * Mz);
  Vector front_layer = createVector(Mx * Mz);
  Vector back_layer = createVector(Mx * Mz);
  Vector top_layer = createVector(Mx * My);
  Vector bottom_layer = createVector(Mx * My);

  // create buffer of the entire room
  Vector AA = NULL;
  if (rank == 0) {
    AA = createVector(Nx * Ny * Nz);
  }
  MPI_Gatherv(A, Mx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, cubes);

  if (rank == 0) {
    printf("Initial:\n");
    printTemperature(AA, Nx, Ny, Nz);
    printf("\n");
  }

  // exchange ghost cells for the very first iteration
  MPI_Isend(A, 1, x_slice, left_rank, 0, cubes, &leftSrequest);
  MPI_Isend(&(A[IDX_3D(Mx - 1, 0, 0, Mx, My)]), 1, x_slice, right_rank, 0, cubes, &rightSrequest);
  MPI_Isend(A, 1, y_slice, front_rank, 0, cubes, &frontSrequest);
  MPI_Isend(&(A[IDX_3D(0, My - 1, 0, Mx, My)]), 1, y_slice, back_rank, 0, cubes, &backSrequest);
  MPI_Isend(A, 1, z_slice, top_rank, 0, cubes, &topSrequest);
  MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Mx, My)]), 1, z_slice, bottom_rank, 0, cubes, &bottomSrequest);

  MPI_Irecv(left_layer, My * Mz, MPI_FLOAT, left_rank, 0, cubes, &leftRrequest);
  MPI_Irecv(right_layer, My * Mz, MPI_FLOAT, right_rank, 0, cubes, &rightRrequest);
  MPI_Irecv(front_layer, Mx * Mz, MPI_FLOAT, front_rank, 0, cubes, &frontRrequest);
  MPI_Irecv(back_layer, Mx * Mz, MPI_FLOAT, back_rank, 0, cubes, &backRrequest);
  MPI_Irecv(top_layer, Mx * My, MPI_FLOAT, top_rank, 0, cubes, &topRrequest);
  MPI_Irecv(bottom_layer, Mx * My, MPI_FLOAT, bottom_rank, 0, cubes, &bottomRrequest);

  // for each time step ..
  for (int t = 0; t < T; t++) {
#ifdef DEBUG
    printf("rank: %d, timestep: %d/%d\n", rank, t, T);
#endif
    // .. we propagate the temperature
    for (int z = 0; z < Mz; z++) {
      // wait for ghost cell exchanges from previous iteration to finish
      // wait for data from neighbouring ranks + wait for sent data to be buffered, because we will change it henceforth
      if (z == 0) {
        MPI_Wait(&topRrequest, MPI_STATUS_IGNORE);
        MPI_Wait(&topSrequest, MPI_STATUS_IGNORE);
      } if (z == Mz - 1) {
        MPI_Wait(&bottomRrequest, MPI_STATUS_IGNORE);
        MPI_Wait(&bottomSrequest, MPI_STATUS_IGNORE);
      }
      for (int y = 0; y < My; y++) {
        if (y == 0 && z == 0) {
          MPI_Wait(&frontRrequest, MPI_STATUS_IGNORE);
          MPI_Wait(&frontSrequest, MPI_STATUS_IGNORE);
        } if (y == My - 1 && z == 0) {
          MPI_Wait(&backRrequest, MPI_STATUS_IGNORE);
          MPI_Wait(&backSrequest, MPI_STATUS_IGNORE);
        }
        for (int x = 0; x < Mx; x++) {
          if (x == 0 && y == 0 && z == 0) {
            MPI_Wait(&leftRrequest, MPI_STATUS_IGNORE);
            MPI_Wait(&leftSrequest, MPI_STATUS_IGNORE);
          } if (x == Mx - 1 && y == 0 && z == 0) {
            MPI_Wait(&rightRrequest, MPI_STATUS_IGNORE);
            MPI_Wait(&rightSrequest, MPI_STATUS_IGNORE);
          }
          // get the current local idx
          int i = IDX_3D(x, y, z, Mx, My);

          // center stays constant (the heat is still on)
          if (local2global(rank, i, Mx, My, Mz, Cx, Cy) == global_source) {
            B[i] = A[i];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i];

          // get temperatures of adjacent cells
          value_t tl = (x != 0) ? A[IDX_3D(x - 1, y, z, Mx, My)] : left_layer[IDX_2D(y, z, My)];
          value_t tr = (x != Mx - 1) ? A[IDX_3D(x + 1, y, z, Mx, My)] : right_layer[IDX_2D(y, z, My)];
          value_t tf = (y != 0) ? A[IDX_3D(x, y - 1, z, Mx, My)] : front_layer[IDX_2D(x, z, Mx)];
          value_t tb = (y != My - 1) ? A[IDX_3D(x, y + 1, z, Mx, My)] : back_layer[IDX_2D(x, z, Mx)];
          value_t tu = (z != 0) ? A[IDX_3D(x, y, z - 1, Mx, My)] : top_layer[IDX_2D(x, y, Mx)];
          value_t td = (z != Mz - 1) ? A[IDX_3D(x, y, z + 1, Mx, My)] : bottom_layer[IDX_2D(x, y, Mx)];

          // compute new temperature at current position
          B[i] = tc + 0.16666 * (tl + tr + tu + td + tf + tb + (-6 * tc));

          // send left/right layer to corresponding ranks immediately after the entire corresponding x_layer has been computed
          // also start receiving from neighbouring ranks - don't block, just continue calculations
          // don't send/receive in last timestep or you will get errors during gather!
          if (x == 0 && y == My -1 && z == Mz - 1 && t != T - 1) {
            MPI_Isend(A, 1, x_slice, left_rank, 0, cubes, &leftSrequest);
            MPI_Irecv(left_layer, My * Mz, MPI_FLOAT, left_rank, 0, cubes, &leftRrequest);
          } if (x == Mx - 1 && y == My - 1 && z == Mz - 1 && t != T - 1) {
            MPI_Isend(&(A[IDX_3D(Mx - 1, 0, 0, Mx, My)]), 1, x_slice, right_rank, 0, cubes, &rightSrequest);
            MPI_Irecv(right_layer, My * Mz, MPI_FLOAT, right_rank, 0, cubes, &rightRrequest);
          }
        }
        if (y == 0 && z == Mz - 1 && t != T - 1) {
          MPI_Isend(A, 1, y_slice, front_rank, 0, cubes, &frontSrequest);
          MPI_Irecv(front_layer, Mx * Mz, MPI_FLOAT, front_rank, 0, cubes, &frontRrequest);
        } if (y == My - 1 && z == Mz - 1 && t != T - 1) {
          MPI_Isend(&(A[IDX_3D(0, My - 1, 1, Mx, My)]), 1, y_slice, back_rank, 0, cubes, &backSrequest);
          MPI_Irecv(back_layer, Mx * Mz, MPI_FLOAT, back_rank, 0, cubes, &backRrequest);
        }
      }
      if (z == 0 && t != T - 1) {
        MPI_Isend(A, 1, z_slice, top_rank, 0, cubes, &topSrequest);
        MPI_Irecv(top_layer, Mx * My, MPI_FLOAT, top_rank, 0, cubes, &topRrequest);
      } if (z == Mz - 1 && t != T - 1) {
        MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Mx, My)]), 1, z_slice, bottom_rank, 0, cubes, &bottomSrequest);
        MPI_Irecv(bottom_layer, Mx * My, MPI_FLOAT, bottom_rank, 0, cubes, &bottomRrequest);
      }
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

#ifdef VERBOSE
    // show intermediate step
    if (0 && !(t % 10000)) {
      MPI_Gatherv(A, Mx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, cubes);

      if (rank == 0) {
        printf("t: %d\n", t);
        printTemperature(AA, Nx, Ny, Nz);
        printf("\n");
      }
    }
#endif
  }

  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);
  releaseVector(front_layer);
  releaseVector(back_layer);
  releaseVector(left_layer);
  releaseVector(right_layer);

#ifdef DEBUG
  printf("hello from rank %d", rank);
#endif
  MPI_Gatherv(A, Mx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, cubes);

  free(recv_count_array);
  free(displacement_array);
  releaseVector(A);

  if (rank == 0) {
#ifdef VERBOSE
    printf("Final:\n");
    printTemperature(AA, Nx, Ny, Nz);
    printf("\n");
#endif

    // ---------- check ----------
    double residual = is_verified_3D(AA, Nx, Ny, Nz, source_x, source_y, source_z, T);
    printf("The maximal deviation from the 1D theory is %fK.", residual);

    // ---------- cleanup ----------
    releaseVector(AA);

    double end = MPI_Wtime();
    printf("The process took %f seconds to finish. \n", end - start);
  }

  MPI_Finalize();

  // done
  return EXIT_SUCCESS;
}
