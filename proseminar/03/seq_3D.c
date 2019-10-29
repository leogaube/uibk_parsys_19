#include <stdio.h>
#include <stdlib.h>

#include <time.h>

typedef float value_t;

#define RESOLUTION 25

// -- vector utilities --

typedef value_t *Vector;

Vector createVector(int N);

void releaseVector(Vector m);

void printTemperature(Vector m, int Nx, int Ny, int Nz);

int MAX(int x, int y){
  return (x > y) ? x : y;
}
int MIN(int x, int y)
{
  return (x < y) ? x : y;
}

int IDX_3D(int x, int y, int z, int nx, int ny)
{
  return x + (y * nx) + (z * nx * ny);
}

// -- simulation code ---

int main(int argc, char **argv) {
  clock_t start = clock();

  // 'parsing' optional input parameter = problem size
  int Nx = 100;
  int Ny = 100;
  int Nz = 100;

  /*if (argc > 3) {
    Nx = atoi(argv[1]);
    Ny = atoi(argv[1]);
    Nz = atoi(argv[1]);
  }*/
  int T = Nx * 500;
  printf("Computing heat-distribution for room size N=%dx%dx%d for T=%d timesteps\n", Nx, Ny, Nz, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx*Ny*Nz);

  // set up initial conditions in A
  for (int i = 0; i < Nx*Ny*Nz; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;
  int source_i = IDX_3D(source_x, source_y, source_z, Nx, Ny);
  A[source_i] = 273 + 60;

  printf("Initial:\n");
  printTemperature(A, Nx, Ny, Nz);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx*Ny*Nz);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (int zi = 0; zi < Nz; zi++)
    {
      for (int yi = 0; yi < Ny; yi++)
      {
        for (int xi = 0; xi < Nx; xi++)
        {
          int i = IDX_3D(xi, yi, zi, Nx, Ny);
          // center stays constant (the heat is still on)
          if (source_i == i)
          {
            B[i] = A[i];
            continue;
          }

          // get temperature at current position
          value_t tc = A[i];

          // get temperatures of adjacent cells
          value_t tl = (xi != 0) ? A[i-1] : tc;
          value_t tr = (xi != Nx - 1) ? A[i+1] : tc;

          value_t tu = (yi != 0) ? A[i-Nx] : tc;
          value_t td = (yi != Ny - 1) ? A[i+Nx] : tc;

          value_t tf = (zi != 0) ? A[i-(Nx*Ny)] : tc;
          value_t tb = (zi != Nz - 1) ? A[i+(Nx * Ny)] : tc;

          // compute new temperature at current position
          B[i] = tc + 0.1666 * (tl + tr + tu + td + tf + tb + (-6 * tc));
        }
      }
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 1000))
    {
      printf("Step t=%d:\n", t);
      printTemperature(A, Nx, Ny, Nz);
      printf("\n");
    }
  }

  releaseVector(B);

  // ---------- check ----------

  printf("Final:\n");
  printTemperature(A, Nx, Ny, Nz);
  printf("\n");

  int success = 1;
  for (int i = 0; i < Nx*Ny; i++) {
    value_t temp = A[i];
    if (273 <= temp && temp <= 273 + 60)
      continue;
    success = 0;
    break;
  }

  printf("Verification: %s\n", (success) ? "OK" : "FAILED");

  // ---------- cleanup ----------

  releaseVector(A);

  clock_t end = clock();
  printf("The process took %g seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);

  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int Nx, int Ny, int Nz) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int Wx = MIN(RESOLUTION, Nx);
  int Wy = MIN(RESOLUTION, Ny);
  int Wz = MIN(RESOLUTION, Nz);

  // step size in each dimension
  int sWx = Nx / Wx;
  int sWy = Ny / Wy;
  int sWz = Nz / Wz;

  // room
  // actual room
  for (int zi = 0; zi < Wz; zi++)
  {
    if (Wz != 1)
      printf("########### SLICE %d ##############\n", zi);
    for (int yi = 0; yi < Wy; yi++)
    {
      // left wall
      printf("X");
      for (int xi = 0; xi < Wx; xi++)
      {
        // get max temperature in this tile
        value_t max_t = 0;
        for (int z = sWz * zi; z < sWz * zi + sWz; z++)
        {
          for (int y = sWy * yi; y < sWy * yi + sWy; y++)
          {
            for (int x = sWx * xi; x < sWx * xi + sWx; x++)
            {
              max_t = (max_t < m[IDX_3D(x, y, z, Nx, Ny)]) ? m[IDX_3D(x, y, z, Nx, Ny)] : max_t;
            }
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
  printf("\n##################################\n");
}
