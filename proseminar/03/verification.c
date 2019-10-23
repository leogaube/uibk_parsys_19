#include <stdlib.h>
#include <stdio.h>

typedef double value_t;
typedef value_t *Vector;

Vector get_result_1D(int source_idx, int room_size);
Vector createVector(int N);
void releaseVector(Vector m);


int main(int argc, char **argv) {
	// TODO get the multidimensional data
	// TODO get how many dims there are
	// TODO get the indices of the source
	int source_idx = 0; // TODO use the first of the source indices
	// TODO get the room size
	int room_size = 100;
	// TODO get the uncertainty below which a deviation is still valid
	double uncert_T = 1e-3;

	Vector result_1D = get_result_1D(source_idx, room_size);

	//TODO go through matrix and compare to result_1D
	/*
	 * The compare goes as follows
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
