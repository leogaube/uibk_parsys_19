#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "heat_stencil.h"

int local2global(int r, int i, int mx, int my, int mz, int cx, int cy){
    int global_x = (r % cx)*mx + (i % mx);
    int global_y = ((r % (cx * cy)) / cx)*my + ((i % (mx * my)) / mx);
    int global_z = (r / (cx * cy)) * mz + (i / (mx * my));

    return IDX_3D(global_x, global_y, global_z, mx * cx, my * cy);
}

Vector createVector(int N){
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int nx, int ny, int nz)
{
#ifdef SLICE_3D
  if (nz != 1){
    nz=1;
    printf("########### ONLY DISPLAYING SLICE 0 of 3D room ##############\n");
  }
#endif

  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int Wx = MIN(RESOLUTION, nx);
  int Wy = MIN(RESOLUTION, ny);
  int Wz = MIN(RESOLUTION, nz);

  // step size in each dimension
  int sWx = nx / Wx;
  int sWy = ny / Wy;
  int sWz = nz / Wz;

  // room
  // actual room
  for (int k = 0; k < Wz; k++)
  {
    if (Wz != 1)
      printf("########### SLICE %d ##############\n", k);
    for (int j = 0; j < Wy; j++)
    {
      // left wall
      printf("X");
      for (int i = 0; i < Wx; i++)
      {
        // get max temperature in this tile
        value_t max_t = 0;
        for (int z = sWz * k; z < sWz * k + sWz; z++)
        {
          for (int y = sWy * j; y < sWy * j + sWy; y++)
          {
            for (int x = sWx * i; x < sWx * i + sWx; x++)
            {
              max_t = (max_t < m[IDX_3D(x, y, z, nx, ny)]) ? m[IDX_3D(x, y, z, nx, ny)] : max_t;
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
  }
  printf("\n############# END OF TIMESTEP ##############\n");
}


double is_verified_3D(Vector result_3D, int nx, int ny, int nz, int source_x, int source_y, int source_z, int T){
	/**
	 * select longest possible distance to the source as 1D room size
	 * as there are no diagonal paths, the norm that is used for this
	 * has to be d((x1,y1),(x2,y2))=|x1-x2|+|y1-y2|
	 * if the euclidian norm would be used to calculate the distance
	 * use:
	 * room_size_1D = (int)ceil(
	 *		sqrt( max_distance_x*max_distance_x +
	 *			  max_distance_y*max_distance_y ) );
	 */
	int max_distance_x = MAX(nx-source_x, source_x);
	int max_distance_y = MAX(ny-source_y, source_y);
	int max_distance_z = MAX(nz-source_z, source_z);
	int room_size_1D = max_distance_x+max_distance_y + max_distance_z;
	Vector result_1D = get_result_1D(room_size_1D, T);

	// Go through the matrix and compare the cells to the result_1D.
	double max_residual = 0;
	for(int z=0; z<nz; z++){
		for(int y=0; y<ny; y++){
			for(int x=0; x<nx; x++){
				double expected_T = get_expected_T(result_1D,
						abs(x-source_x), abs(y-source_y), abs(z-source_z), room_size_1D);

				double residual = result_3D[IDX_3D(x,y,z,nx,ny)]-expected_T;
				// The residual can be compared to a maximally allowed value here.
				// As this maximal value depends on the chosen geometry/norm, the
				// maximal residual is found and returned instead.
				if(abs(residual)>abs(max_residual)){
					max_residual = residual;
				}
			}
		}
	}
	return max_residual;
}


double is_verified_2D(Vector result_2D, int nx, int ny, int source_x, int source_y, int T){
	return is_verified_3D(result_2D, nx, ny, 1, source_x, source_y, 0, T);
}


/**
 * This function gets a 1D result with the code from the prev. homework.
 * The source is on idx 0.
 */
Vector get_result_1D(int N, int T){
  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = createVector(N);

  // set up initial conditions in A
  for (int i = 0; i < N; i++) {
	A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source at idx 0
  int source_x = 0;
  A[source_x] = 273 + 60;
#ifdef VERBOSE
	  printf("1D initial:\n");
	  printTemperature(A, N, 1, 1);
#endif

  // ---------- compute ----------
  // create a second buffer for the computation
  Vector B = createVector(N);

  // for each time step ..
  for (int t = 0; t < T; t++) {
	// .. we propagate the temperature
	for (long long i = 0; i < N; i++) {
	  // center stays constant (the heat is still on)
	  if (i == source_x) {
		B[i] = A[i];
		continue;
	  }

	  // get temperature at current position
	  value_t tc = A[i];

	  // get temperatures of adjacent cells
	  value_t tl = (i != 0) ? A[i - 1] : tc;
	  value_t tr = (i != N - 1) ? A[i + 1] : tc;

	  // compute new temperature at current position
	  B[i] = tc + 0.2 * (tl + tr + (-2 * tc));
	}

	// swap matrices (just pointers, not content)
	Vector H = A;
	A = B;
	B = H;
  }

  releaseVector(B);
#ifdef VERBOSE
	  printf("1D final:\n");
	  printTemperature(A, N, 1, 1);
#endif

  return A;
}


/*
 * The compare goes as follows:
 * 1. calculate the distance d to the source from the indices
 *    (The norm to calculate this distance is not the euclidian norm as the
 *    distance there is no diagonal path between cells.)
 * 2. round the distance up and down to get the closest cells d_1 and d_2
 *    from the result_1D vector
 * 3. the temperature expected from the result_1D is an interpolation
 *    between the closest cells T_1 and T_2, using
 *    T = T_1 + (T_2 - T_1) * (d_2 - d)
 *    where the index 1 is the cell closer to the source
 * 4. The difference between T and the N-dim cell's temperature has to
 *    be below the uncert_T.
 */
double get_expected_T(Vector result_1D, int dx, int dy, int dz, int room_size_1D){
	double d = dx+dy+dz; // euclidian would be sqrt(dx*dx+dy*dy+dz*dz);
	double T_1 = result_1D[(int)floor(d)];
	double T_2 = result_1D[(int)ceil(d)];
	double T_expected = T_1 + (T_2 - T_1) * (ceil(d) - d);
	return T_expected;
}

