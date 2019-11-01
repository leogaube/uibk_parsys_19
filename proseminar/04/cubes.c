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
  int periods[1] = {0, 0, 0};

  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &cubes);
  MPI_Comm_rank(cubes, &rank);

#ifdef VERBOSE
  if (rank == 0)
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
#endif

  // get the adjacent slices
  int left_rank = rank;
  int right_rank = rank;
  int front_rank = rank;
  int back_rank = rank;
  int top_rank = rank;
  int bottom_rank = rank;

  MPI_Cart_shift(cubes, 0, 1, &left_rank, &right_rank);
  MPI_Cart_shift(cubes, 1, 1, &front_rank, &back_rank);
  MPI_Cart_shift(cubes, 2, 1, &bottom_rank, &top_rank);

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

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Mx * My * Mz);

  // set up initial conditions in A
  for (int i = 0; i < Mx * My * Mz; i++) {
    A[i] = 273;  // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;
  int source_i = IDX_3D(source_x, source_y, source_z, Nx, Ny);

  if (rank * Mx <= source_x && (rank + 1) * Mx > source_x &&
      rank * My <= source_y && (rank + 1) * My > source_y &&
      rank * Mz <= source_z && (rank + 1) * Mz > source_z) {
    printf("hello from source init");
    A[IDX_3D(source_x - (rank * Mx), source_y - (rank * My), source_z - (rank * Mz), Nx, Ny)] = 273 + 60;
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Mx * My * Mz);

  // create ghost exchange layers
  Vector top_layer = createVector(Nx * Ny);
  Vector bottom_layer = createVector(Nx * Ny);
  Vector front_layer = createVector(Nx * Nz);
  Vector back_layer = createVector(Nx * Nz);
  Vector left_layer = createVector(Ny * Nz);
  Vector right_layer = createVector(Ny * Nz);

  // create buffer of the entire room
  Vector AA = NULL;
  if (rank == 0) {
    AA = createVector(Nx * Ny * Nz);
  }

  // exchange ghost cells for the very first iteration
  MPI_Isend(A, Nx * Ny, MPI_FLOAT, top_rank, 0, cubes, &topSrequest);
  MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Nx, Ny)]), Nx * Ny, MPI_FLOAT, bottom_rank, 0, cubes, &bottomSrequest);
  MPI_Isend(A, Nx * Ny, MPI_FLOAT, top_rank, 0, cubes, &topSrequest);
  MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Nx, Ny)]), Nx * Ny, MPI_FLOAT, bottom_rank, 0, cubes, &bottomSrequest);
  MPI_Isend(A, Nx * Ny, MPI_FLOAT, top_rank, 0, cubes, &topSrequest);
  MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Nx, Ny)]), Nx * Ny, MPI_FLOAT, bottom_rank, 0, cubes, &bottomSrequest);

  MPI_Irecv(top_layer, Nx * Ny, MPI_FLOAT, top_rank, 0, cubes, &topRrequest);
  MPI_Irecv(bottom_layer, Nx * Ny, MPI_FLOAT, bottom_rank, 0, cubes, &bottomRrequest);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (int z = 0; z < Mz; z++) {
      // wait for ghost cell exchanges from previous iteration to finish
      // wait for data from neighbouring ranks + wait for sent data to be buffered, because we will change it henceforth
      if (z == 0) {
        MPI_Wait(&topRrequest, MPI_STATUS_IGNORE);
        MPI_Wait(&topSrequest, MPI_STATUS_IGNORE);
      } else if (z == Mz - 1) {
        MPI_Wait(&bottomRrequest, MPI_STATUS_IGNORE);
        MPI_Wait(&bottomSrequest, MPI_STATUS_IGNORE);
      }
      for (int y = 0; y < My; y++) {
        if (y == 0) {
          MPI_Wait(&frontRrequest, MPI_STATUS_IGNORE);
          MPI_Wait(&frontSrequest, MPI_STATUS_IGNORE);
        } else if (y == My - 1) {
          MPI_Wait(&backRrequest, MPI_STATUS_IGNORE);
          MPI_Wait(&backSrequest, MPI_STATUS_IGNORE);
        }
        for (int x = 0; x < Mx; x++) {
          if (y == 0) {
            MPI_Wait(&leftRrequest, MPI_STATUS_IGNORE);
            MPI_Wait(&leftSrequest, MPI_STATUS_IGNORE);
          } else if (y == My - 1) {
            MPI_Wait(&rightRrequest, MPI_STATUS_IGNORE);
            MPI_Wait(&rightSrequest, MPI_STATUS_IGNORE);
          }
          // get the current local idx
          int i = IDX_3D(x, y, z, Mx, My);

          // center stays constant (the heat is still on)
          if (local2global(rank, i, Mx, My, Mz, Cx, Cy) == IDX_3D(source_x, source_y, source_z, Nx, Ny)) {
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

          //printf("%f, %f, %f, %f, %f, %f, %f", tc, tl, tr, tu, td, tf, tb);

          // compute new temperature at current position
          B[i] = tc + 0.16666 * (tl + tr + tu + td + tf + tb + (-6 * tc));
        }
      }
      // send the uppermost/lowest layer to upper/lower rank right after the entire z-layer has been computed
      // also start receiving from neighbouring ranks - don't block, just continue calculations
      // don't send/receive in last timestep or you will get errors during gather!
      if (z == 0 && t != T - 1) {
        MPI_Isend(A, Nx * Ny, MPI_FLOAT, top_rank, 0, cubes, &topSrequest);
        MPI_Irecv(top_layer, Nx * Ny, MPI_FLOAT, top_rank, 0, cubes, &topRrequest);
      } else if (z == Mz - 1 && t != T - 1) {
        MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Nx, Ny)]), Nx * Ny, MPI_FLOAT, bottom_rank, 0, cubes, &bottomSrequest);
        MPI_Irecv(bottom_layer, Nx * Ny, MPI_FLOAT, bottom_rank, 0, cubes, &bottomRrequest);
      }
    }
    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

#ifdef VERBOSE
    // show intermediate step
    if (!(t % 1000)) {
      MPI_Gather(A, Nx * Ny * Mz, MPI_FLOAT, AA, Nx * Ny * Mz, MPI_FLOAT, 0, cubes);
      if (rank == 0) {
        printf("Step t=%d:\n", t);
        printTemperature(AA, Nx, Ny, Nz);
        printf("\n");
      }
    }
#endif
  }

  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);

  MPI_Gather(A, Nx * Ny * Mz, MPI_FLOAT, AA, Nx * Ny * Mz, MPI_FLOAT, 0, cubes);

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
