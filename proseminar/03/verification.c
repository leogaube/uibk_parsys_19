#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "heat_stencil.h"

Vector get_result_1D(int N, int T);
double get_expected_T(Vector result_1D, int x, int y, int room_size_1D);

int is_verified_2D(Vector result_2D, int nx, int ny, int source_x, int source_y, int T){
	double uncert_T = 1e-3;

	// select longest possible distance to the source as 1D room size
	int max_distance_x = MAX(nx-source_x, source_x);
	int max_distance_y = MAX(ny-source_y, source_y);
	int room_size_1D = (int)ceil(
			sqrt( max_distance_x*max_distance_x +
				  max_distance_y*max_distance_y ) );
	Vector result_1D = get_result_1D(room_size_1D, T);

	// Go through the matrix and compare the cells to the result_1D.
	for(int y=0; y<ny; y++){
		for(int x=0; x<nx; x++){
			double expected_T = get_expected_T(result_1D,
					x-source_x, y-source_y, room_size_1D);

			double residual = abs(result_2D[IDX_2D(x,y,nx)]-expected_T);

			if(residual>=uncert_T){
				printf("The difference is too large (%f K) for the "
						"cell with the indices (%i,%i).\n"
						"Expected: %f K, Real: %f K.\n",
						residual, x, y, expected_T, result_2D[IDX_2D(x,y,nx)]);
				return -1;
			}
		}
	}
	return 0;
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
	  return A;
}


/*
 * The compare goes as follows:
 * 1. calculate the euclidian distance d to the source from the indices
 * 2. round the distance up and down to get the closest cells d_1 and d_2
 *    from the result_1D vector
 * 3. the temperature expected from the result_1D is an interpolation
 *    between the closest cells T_1 and T_2, using
 *    T = T_1 + (T_2 - T_1) * (d_2 - d)
 *    where the index 1 is the cell closer to the source
 * 4. The difference between T and the N-dim cell's temperature has to
 *    be below the uncert_T.
 */
double get_expected_T(Vector result_1D, int x, int y, int room_size_1D){
	double d = sqrt(x*x+y*y);
	double T_1 = result_1D[(int)floor(d)];
	double T_2 = result_1D[(int)ceil(d)];
	double T_expected = T_1 + (T_2 - T_1) * (ceil(d) - d);
	return T_expected;
}
