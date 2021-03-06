#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "heat_stencil.h"

// -- simulation code ---

int main(int argc, char **argv)
{
  double start = MPI_Wtime();

  // 'parsing' optional input parameter = problem size
  int Nx = 32;
  int Ny = 32;
  int Nz = 32;
  if (argc == 2)
  {
    Nx = Ny = Nz = atoi(argv[1]);
  }
  else if (argc == 4)
  {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    Nz = atoi(argv[3]);
  }
  int T = MAX(MAX(Nx, Ny), Nz) * TIMESTEPS_MUL;

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
  int Mz = Nz / numProcs;

  MPI_Request topRrequest;
  MPI_Request topSrequest;
  MPI_Request bottomRrequest;
  MPI_Request bottomSrequest;
  MPI_Comm slabs;
  int dims[1] = {numProcs};
  int periods[1] = {0};
  MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &slabs);
  MPI_Comm_rank(slabs, &rank);

#ifdef VERBOSE
  if (rank == 0)
    printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
#endif

  // get the adjacent slices
  int top_rank = rank;
  int bottom_rank = rank;
  MPI_Cart_shift(slabs, 0, 1, &bottom_rank, &top_rank);
  if (top_rank == MPI_PROC_NULL)
  {
    top_rank = rank;
  }
  if (bottom_rank == MPI_PROC_NULL)
  {
    bottom_rank = rank;
  }

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx * Ny * Mz);

  // set up initial conditions in A
  for (long i = 0; i < Nx * Ny * Mz; i++)
  {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;

  if (rank * Mz <= source_z && (rank + 1) * Mz > source_z)
  {
    A[IDX_3D(source_x, source_y, source_z - (rank * Mz), Nx, Ny)] = 273 + 60;
  }

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx * Ny * Mz);
  Vector top_layer = createVector(Nx * Ny);
  Vector bottom_layer = createVector(Nx * Ny);

  // create buffer of the entire room
  Vector AA = NULL;
  if (rank == 0)
  {
    AA = createVector(Nx * Ny * Nz);
  }

  // exchange ghost cells for the very first iteration
  MPI_Isend(A, Nx * Ny, MPI_FLOAT, top_rank, 0, slabs, &topSrequest);
  MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Nx, Ny)]), Nx * Ny, MPI_FLOAT, bottom_rank, 0, slabs, &bottomSrequest);

  MPI_Irecv(top_layer, Nx * Ny, MPI_FLOAT, top_rank, 0, slabs, &topRrequest);
  MPI_Irecv(bottom_layer, Nx * Ny, MPI_FLOAT, bottom_rank, 0, slabs, &bottomRrequest);

  // for each time step ..
  for (int t = 0; t < T; t++)
  {
    // .. we propagate the temperature
    for (int z = 0; z < Mz; z++)
    {
      // wait for ghost cell exchanges from previous iteration to finish
      // wait for data from neighbouring ranks + wait for sent data to be buffered, because we will change it henceforth
      if (z == 0)
      {
        MPI_Wait(&topRrequest, MPI_STATUS_IGNORE);
        MPI_Wait(&topSrequest, MPI_STATUS_IGNORE);
      }
      if (z == Mz - 1)
      {
        MPI_Wait(&bottomRrequest, MPI_STATUS_IGNORE);
        MPI_Wait(&bottomSrequest, MPI_STATUS_IGNORE);
      }
      for (int y = 0; y < Ny; y++)
      {
        for (int x = 0; x < Nx; x++)
        {
          // get the current local idx
          long i = IDX_3D(x, y, z, Nx, Ny);

          // center stays constant (the heat is still on)
          if (IDX_3D(x, y, z + (rank*Mz), Nx, Ny) == IDX_3D(source_x, source_y, source_z, Nx, Ny))
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

          //printf("%f, %f, %f, %f, %f, %f, %f", tc, tl, tr, tu, td, tf, tb);

          // compute new temperature at current position
          B[i] = tc + 0.16666 * (tl + tr + tu + td + tf + tb + (-6 * tc));
        }
      }
      // send the uppermost/lowest layer to upper/lower rank right after the entire z-layer has been computed
      // also start receiving from neighbouring ranks - don't block, just continue calculations
      // don't send/receive in last timestep or you will get errors during gather!
      if (z == 0 && t != T - 1)
      {
        MPI_Isend(A, Nx * Ny, MPI_FLOAT, top_rank, 0, slabs, &topSrequest);
        MPI_Irecv(top_layer, Nx * Ny, MPI_FLOAT, top_rank, 0, slabs, &topRrequest);
      }
      if (z == Mz - 1 && t != T-1)
      {
        MPI_Isend(&(A[IDX_3D(0, 0, Mz - 1, Nx, Ny)]), Nx * Ny, MPI_FLOAT, bottom_rank, 0, slabs, &bottomSrequest);
        MPI_Irecv(bottom_layer, Nx * Ny, MPI_FLOAT, bottom_rank, 0, slabs, &bottomRrequest);
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
      MPI_Gather(A, Nx * Ny * Mz, MPI_FLOAT, AA, Nx * Ny * Mz, MPI_FLOAT, 0, slabs);
      if (rank == 0)
      {
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

  MPI_Gather(A, Nx * Ny * Mz, MPI_FLOAT, AA, Nx * Ny * Mz, MPI_FLOAT, 0, slabs);
  
  releaseVector(A);

  if (rank == 0)
  {
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
