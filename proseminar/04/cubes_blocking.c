#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "heat_stencil.h"

// -- simulation code ---
int main(int argc, char** argv) {
  double start = MPI_Wtime();

  // 'parsing' optional input parameters: room size (&& rank layout)
  int Nx = 32;
  int Ny = 32;
  int Nz = 32;
  int Px, Py, Pz;
  int rank_layout = 0;
  if (argc == 2) {
    Nx = Ny = Nz = atoi(argv[1]);
  } else if (argc == 3) {
    Nx = Ny = Nz = atoi(argv[1]);
    Px = Py = Pz = atoi(argv[2]);
  } else if (argc == 4) {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
  } else if (argc == 7) {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
    Px = atoi(argv[4]);
    Py = atoi(argv[5]);
    Pz = atoi(argv[6]);
    rank_layout = 1;
  } else {
    printf("wrong number of arguments!");
    return EXIT_FAILURE;
  }

  if (Nx < Ny || Ny < Nz) {
    printf("Please specify a room that meets the contition: Nx >= Ny >= Nz!");
    return EXIT_FAILURE;
  }
  int T = MAX(MAX(Nx, Ny), Nz) * TIMESTEPS_MUL;

  // MPI setup
  int rank, numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (Nx * Ny * Nz < numProcs) {
    printf("Room is too small for %d ranks", numProcs);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  if (rank_layout) {
    if (numProcs != Px * Py * Pz) {
      printf("#ranks does not match with provided Px, Py and Pz!");
      MPI_Finalize();
      return EXIT_FAILURE;
    }
  } else {
    // find the optimal rank layout
    // e.g. 8 slots --> Px=2, Py=2 Pz=2
    // e.g. 42 slots --> Px=3, Py=2 Pz=7
    // e.g. 17 slots --> Px=1, Py=1, Pz=17
    int c = floorf(cbrt(numProcs));
    while (numProcs % c)
      c--;
    Px = c;

    int s = floorf(sqrt(numProcs / Px));
    while ((numProcs / Px) % s)
      s--;
    Py = s;
    Pz = numProcs / (Px * Py);

    // also allow non-cubic 1D, 2D or 3D room sizes
    if (Nz % Pz) {
      Px *= Pz;
      Pz = 1;
    }
    if (Ny % Py) {
      Px *= Py;
      Py = 1;
    }
  }

  // same program works for poles and cubes --> identify by name of executable
  if (strstr(argv[0], "poles_blocking") != NULL) {
    if (Px != 1 && Py != 1 && Pz != 1) {
      Px *= Py;
      Py = 1;
    }
  }

  if (Nx % Px || Ny % Py || Nz % Pz) {
    printf("Room size (%dx%dx%d) cannot be split up evenly among ranks (%dx%dx%d)\n", Nx, Ny, Nz, Px, Py, Pz);
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  // subroom sizes in all dimensions
  int Mx = Nx / Px;
  int My = Ny / Py;
  int Mz = Nz / Pz;

  // create new communicator
  MPI_Comm cubes;

  int dims[3] = {Px, Py, Pz};
  int periods[3] = {0, 0, 0};

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cubes);
  MPI_Comm_rank(cubes, &rank);

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

  // possibly send data to yourself for simplicity
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

  // 2D subroom slices for each dimension --> useful for sending 2D ghost cell slices to neighbouring ranks
  MPI_Datatype x_slice, y_slice, z_slice;
  MPI_Type_vector(My * Mz, 1, Mx, MPI_FLOAT, &x_slice);
  MPI_Type_vector(Mz, Mx, Mx * My, MPI_FLOAT, &y_slice);
  MPI_Type_vector(1, Mx * My, 1, MPI_FLOAT, &z_slice);
  MPI_Type_commit(&x_slice);
  MPI_Type_commit(&y_slice);
  MPI_Type_commit(&z_slice);

  // create subroom_datatypes for MPI_Gatherv
  //--> start of subroom determined by displacement_array, pretend size of subroom datatype == 1*sizeof(value_t)
  int room_sizes[3] = {Nx, Ny, Nz};
  int subroom_sizes[3] = {Mx, My, Mz};  // for some reason Mz and Mx have to be swaped or MPI_ORDER_FORTRAN?!
  int start_array[3] = {0, 0, 0};
  MPI_Datatype cubic_subroom, resized_subroom;
  MPI_Type_create_subarray(3, room_sizes, subroom_sizes, start_array, MPI_ORDER_FORTRAN, MPI_FLOAT, &cubic_subroom);
  // 'pretend' that subroom is only 1 float in size --> displacements ^= local2global conversion
  MPI_Type_create_resized(cubic_subroom, 0, 1 * sizeof(value_t), &resized_subroom);
  MPI_Type_commit(&resized_subroom);

  int* recv_count_array = malloc(sizeof(int) * numProcs);
  int* displacement_array = malloc(sizeof(int) * numProcs);
  for (int r = 0; r < numProcs; r++) {
    recv_count_array[r] = 1;
    // get the global index of every rank at its local position 0
    displacement_array[r] = local2global(r, 0, Mx, My, Mz, Px, Py);
  }

  // ---------- setup ----------

  // create a buffer for storing temperature fields for sub-rooms
  Vector A = createVector(Mx * My * Mz);

  // set up initial conditions in A
  for (int i = 0; i < Mx * My * Mz; i++) {
    A[i] = 273;  // debug with: + ((double)rank/(numProcs+1)) * 60;  // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;
  int global_source = IDX_3D(source_x, source_y, source_z, Nx, Ny);
  int local_source = IDX_3D(source_x % Mx, source_y % My, source_z % Mz, Mx, My);

  if (local2global(rank, local_source, Mx, My, Mz, Px, Py) == global_source) {
    A[local_source] = 273 + 60;
  }

  // ---------- compute ----------

  // create buffer of the entire room
  Vector AA = NULL;
  if (rank == 0) {
    AA = createVector(Nx * Ny * Nz);
  }

#ifdef VERBOSE
  MPI_Gatherv(A, Mx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, cubes);
  if (rank == 0) {
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
    printf("There are %d subrooms (Px=%d, Py=%d, Pz=%d) each of size Mx=%d, My=%d, Mz=%d\n\n", numProcs, Px, Py, Pz, Mx, My, Mz);

    printf("Initial:\n");
    printTemperature(AA, Nx, Ny, Nz);
    printf("\n");
  }
#endif

  // create a second buffer for the computation
  Vector B = createVector(Mx * My * Mz);

  // create ghost exchange layers
  Vector left_layer = createVector(My * Mz);
  Vector right_layer = createVector(My * Mz);
  Vector front_layer = createVector(Mx * Mz);
  Vector back_layer = createVector(Mx * Mz);
  Vector top_layer = createVector(Mx * My);
  Vector bottom_layer = createVector(Mx * My);

  // exchange ghost cells for the very first iteration
  int coords[3];
  MPI_Cart_coords(cubes, rank, 3, coords);
  //printf("rank: %d, cord[0]: %d, cord[1]: %d, cord[2]: %d\n", rank, coords[0], coords[1], coords[2]);
  //printf("rank: %d, l: %d, r: %d, t: %d, b: %d, f: %d, b: %d\n", rank, left_rank, right_rank, top_rank, bottom_rank, front_rank, back_rank);
  // for each time step ..
  for (int t = 0; t < T; t++) {
    if (coords[2] % 2) {
      if (coords[2] != 0) {
        MPI_Ssend(A, 1, x_slice, left_rank, 0, cubes);
        MPI_Recv(left_layer, My * Mz, MPI_FLOAT, left_rank, 0, cubes, MPI_STATUS_IGNORE);
      }
    } else {
      if (coords[2] != Pz - 1) {
        MPI_Recv(right_layer, My * Mz, MPI_FLOAT, right_rank, 0, cubes, MPI_STATUS_IGNORE);
        MPI_Ssend(&(A[IDX_3D(Mx - 1, 0, 0, Mx, My)]), 1, x_slice, right_rank, 0, cubes);
      }
    }
    if (coords[0] % 2) {
      if (coords[0] != 0) {
        MPI_Ssend(A, 1, z_slice, top_rank, 0, cubes);
        MPI_Recv(top_layer, Mx * My, MPI_FLOAT, top_rank, 0, cubes, MPI_STATUS_IGNORE);
      }
    } else {
      if (coords[0] != Px - 1) {
        MPI_Recv(bottom_layer, Mx * My, MPI_FLOAT, bottom_rank, 0, cubes, MPI_STATUS_IGNORE);
        MPI_Ssend(&(A[IDX_3D(0, 0, Mz - 1, Mx, My)]), 1, z_slice, bottom_rank, 0, cubes);
      }
    }
    if (coords[1] % 2) {
      if (coords[1] != 0) {
        MPI_Ssend(A, 1, y_slice, front_rank, 0, cubes);
        MPI_Recv(front_layer, Mx * Mz, MPI_FLOAT, front_rank, 0, cubes, MPI_STATUS_IGNORE);
      }
    } else {
      if (coords[1] != Py - 1) {
        MPI_Recv(back_layer, Mx * Mz, MPI_FLOAT, back_rank, 0, cubes, MPI_STATUS_IGNORE);
        MPI_Ssend(&(A[IDX_3D(0, My - 1, 0, Mx, My)]), 1, y_slice, back_rank, 0, cubes);
      }
    }

    // .. we propagate the temperature
    for (int z = 0; z < Mz; z++) {
      for (int y = 0; y < My; y++) {
        for (int x = 0; x < Mx; x++) {
          // get the current local idx
          int i = IDX_3D(x, y, z, Mx, My);

          // center stays constant (the heat is still on)
          if (local2global(rank, i, Mx, My, Mz, Px, Py) == global_source) {
            B[i] = A[i];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i];

          // get temperatures of adjacent cells
          value_t tl = (rank != left_rank || x != 0) ? (x != 0) ? A[IDX_3D(x - 1, y, z, Mx, My)] : left_layer[IDX_2D(y, z, My)] : tc;
          value_t tr = (rank != right_rank || x != Mx - 1) ? (x != Mx - 1) ? A[IDX_3D(x + 1, y, z, Mx, My)] : right_layer[IDX_2D(y, z, My)] : tc;
          value_t tf = (rank != front_rank || y != 0) ? (y != 0) ? A[IDX_3D(x, y - 1, z, Mx, My)] : front_layer[IDX_2D(x, z, Mx)] : tc;
          value_t tb = (rank != back_rank || y != My - 1) ? (y != My - 1) ? A[IDX_3D(x, y + 1, z, Mx, My)] : back_layer[IDX_2D(x, z, Mx)] : tc;
          value_t tu = (rank != top_rank || z != 0) ? (z != 0) ? A[IDX_3D(x, y, z - 1, Mx, My)] : top_layer[IDX_2D(x, y, Mx)] : tc;
          value_t td = (rank != bottom_rank || z != Mz - 1) ? (z != Mz - 1) ? A[IDX_3D(x, y, z + 1, Mx, My)] : bottom_layer[IDX_2D(x, y, Mx)] : tc;

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
    if (!(t % 1000)) {
      MPI_Gatherv(A, Mx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, cubes);

      if (rank == 0) {
        printf("t: %d\n", t);
        printTemperature(AA, Nx, Ny, Nz);
        printf("\n");
      }
    }
#endif
  }

  MPI_Gatherv(A, Mx * My * Mz, MPI_FLOAT, AA, recv_count_array, displacement_array, resized_subroom, 0, cubes);

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

  //cleanup
  releaseVector(A);
  releaseVector(B);
  releaseVector(bottom_layer);
  releaseVector(top_layer);
  releaseVector(front_layer);
  releaseVector(back_layer);
  releaseVector(left_layer);
  releaseVector(right_layer);

  MPI_Type_free(&x_slice);
  MPI_Type_free(&y_slice);
  MPI_Type_free(&z_slice);
  MPI_Type_free(&cubic_subroom);
  MPI_Type_free(&resized_subroom);

  free(recv_count_array);
  free(displacement_array);

  MPI_Finalize();

  // done
  return EXIT_SUCCESS;
}
