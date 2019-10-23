#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef double value_t;
typedef value_t *Vector;

Vector get_result_1D(int source_idx, int room_size);
Vector createVector(int N);
void releaseVector(Vector m);
double get_expected_T(Vector result_1D, int x, int y, int z, int source_idx, int room_size_1D);

int main(int argc, char **argv) {
	// TODO get the multidimensional data
	double result_3D[10][10][1] = {};
	// TODO get how many dims there are
	int n_dim = 2;
	// TODO get the indices of the source
	int source_x = 1;
	int source_y = 1;
	int source_z = 0; // set zero if 2 dim
	int source_idx = source_x;
	// TODO get the room size
	int room_size = 10;
	// TODO get the uncertainty below which a deviation is still valid
	double uncert_T = 1e-3;

	Vector result_1D = get_result_1D(source_idx, room_size*n_dim);

	// Go through the matrix and compare the cells to the result_1D.
	for(int z=0; z<1+(n_dim-2)*(room_size-1); z++){
		for(int y=0; y<room_size; y++){
			for(int x=0; x<room_size; x++){
				double expected_T = get_expected_T(result_1D,
						x-source_x, y-source_y, z-source_z,
						source_idx, room_size*n_dim);

				double residual = abs(result_3D[x][y][z]-expected_T);

				if(residual>=uncert_T){
					printf("The difference is too large (%f K) for the "
							"cell with the indices (%i,%i,%i).\n"
							"Expected: %f K, Real: %f K.\n",
							residual, x, y, z, expected_T, result_3D[x][y][z]);
					return EXIT_SUCCESS;
				}
			}
		}
	}

	return EXIT_SUCCESS;
}



Vector get_result_1D(int source_idx, int room_size){
	// TODO calculate from the code used in previous homework
	Vector result_1D = createVector(room_size);
	for(int i=0; i<room_size; i++){
		result_1D[i]=273;
	}
	return result_1D;
}

Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

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
double get_expected_T(Vector result_1D, int x, int y, int z, int source_idx, int room_size_1D){
	double d = sqrt(x*x+y*y+z*z);
	double T_1 = result_1D[((int)floor(d)+source_idx)%room_size_1D];
	double T_2 = result_1D[((int)ceil(d)+source_idx)%room_size_1D];
	double T_expected = T_1 + (T_2 - T_1) * (ceil(d) - d);
	return T_expected;
}
