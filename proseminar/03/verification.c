#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "heat_stencil.h"

Vector get_result_1D(int room_size);
double get_expected_T(Vector result_1D, int x, int y, int room_size_1D);

int verify_2D(Vector result_2D, int nx, int ny, int source_x, int source_y){
	double uncert_T = 1e-3;

	// select longest possible distance to source as room size
	int room_size_1D = (int)ceil(
			sqrt( (nx-source_x)*(nx-source_x) +
				  (ny-source_y)*(ny-source_y) ) );
	Vector result_1D = get_result_1D(room_size_1D);

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


Vector get_result_1D(int room_size){
	// TODO calculate from the code used in previous homework
	Vector result_1D = createVector(room_size);
	for(int i=0; i<room_size; i++){
		result_1D[i]=273;
	}
	return result_1D;
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
double get_expected_T(Vector result_1D, int x, int y, int z, int room_size_1D){
	double d = sqrt(x*x+y*y+z*z);
	double T_1 = result_1D[(int)floor(d)];
	double T_2 = result_1D[(int)ceil(d)];
	double T_expected = T_1 + (T_2 - T_1) * (ceil(d) - d);
	return T_expected;
}
