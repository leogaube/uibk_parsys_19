#include <stdio.h>
#include <stdlib.h>

#include <time.h>

#include "heat_stencil.h"

#define RESOLUTION 120

void printTemperature(Vector m, int nx, int ny, int nz);

// -- simulation code ---

int main(int argc, char **argv) {
  clock_t start = clock();

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
  int T = Nx+Ny+Nz;
#ifdef VERBOSE
  printf("Computing heat-distribution for room size Nx=%d, Ny=%d, Nz=%d for T=%d timesteps\n", Nx, Ny, Nz, T);
#endif

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(Nx*Ny*Nz);

  // set up initial conditions in A
  for (long i = 0; i < Nx*Ny*Nz; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = Nx / 4;
  int source_y = Ny / 4;
  int source_z = Nz / 4;
  A[IDX_3D(source_x,source_y, source_z,Nx,Ny)] = 273 + 60;

#ifdef VERBOSE
  printf("Initial:\n");
  printTemperature(A, Nx, Ny, Nz);
#endif
  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = createVector(Nx*Ny*Nz);

  // for each time step ..
  for (int t = 0; t < T; t++) {
	// .. we propagate the temperature
	for (int z = 0; z < Nz; z++) {
		for (int y = 0; y < Ny; y++) {
			for(int x = 0; x < Nx; x++){
				// get the current idx
				long i = IDX_3D(x,y,z,Nx,Ny);

				// center stays constant (the heat is still on)
				if (i == IDX_3D(source_x,source_y,source_z,Nx,Ny)) {
					B[i] = A[i];
					continue;
				}

				// get temperature at current position
				value_t tc = A[i];

				// get temperatures of adjacent cells
				value_t tl = (x != 0) ? A[IDX_3D(x-1,y,z,Nx,Ny)] : tc;
				value_t tr = (x != Nx - 1) ? A[IDX_3D(x+1,y,z,Nx,Ny)] : tc;
				value_t tu = (y != 0) ? A[IDX_3D(x,y-1,z,Nx,Ny)] : tc;
				value_t td = (y != Ny - 1) ? A[IDX_3D(x,y+1,z,Nx,Ny)] : tc;
				value_t tf = (z != 0) ? A[IDX_3D(x,y,z-1,Nx,Ny)] : tc;
				value_t tb = (z != Nz - 1) ? A[IDX_3D(x,y,z+1,Nx,Ny)] : tc;


				// compute new temperature at current position
				B[i] = tc + 0.2 * (tl + tr + tu + td + tf + tb + (-6 * tc));
			}
		}
	}
	// swap matrices (just pointers, not content)
	Vector H = A;
	A = B;
	B = H;
  }

  releaseVector(B);
#ifdef VERBOSE
  printf("Final:\n");
  printTemperature(A, Nx, Ny, Nz);
#endif

  // ---------- check ----------
  double residual = is_verified_3D(A, Nx, Ny, Nz, source_x, source_y, source_z, T);
  printf("The maximal deviation from the 1D theory is %fK.", residual);

  // ---------- cleanup ----------

  releaseVector(A);

  clock_t end = clock();
  printf("The process took %g seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);

  // done
  return EXIT_SUCCESS;
}

void printTemperature(Vector m, int nx, int ny, int nz) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int Wx = RESOLUTION;
  int Wy = RESOLUTION;
  int Wz = RESOLUTION;
  if(Wx>nx){	Wx=nx;	};
  if(Wy>ny){	Wy=ny;	};
  if(Wz>nz){	Wz=nz;	};

  // step size in each dimension
  int sWx = nx / Wx;
  int sWy = ny / Wy;
  int sWz = nz / Wz;

  // room
  // actual room
  for(int k=0; k<Wz; k++){
	  for(int j=0; j<Wy; j++){
		  // left wall
		  printf("X");
		  for (int i = 0; i < Wx; i++) {
			// get max temperature in this tile
			value_t max_t = 0;
			for (int z = sWz * k; z < sWz * k + sWz; z++) {
				for (int y = sWy * j; y < sWy * j + sWy; y++) {
					for (int x = sWx * i; x < sWx * i + sWx; x++) {
						max_t = (max_t < m[IDX_3D(x,y,z,nx,ny)]) ? m[IDX_3D(x,y,z,nx,ny)] : max_t;
					}
				}
			}
			value_t temp = max_t;

			// pick the 'color'
			int c = ((temp - min) / (max - min)) * numColors;
			c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

			// print the average temperature
		    // if numbers are desired use printf("%2.2f\t",temp-273);
			printf("%c", colors[c]);
		  }
		  // right wall
		  printf("X\n");
	  }
	  printf("\n");
  }
}
