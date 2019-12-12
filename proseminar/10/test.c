#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    #pragma omp parallel for
    for(int i = 0; i < omp_get_max_threads(); i++) {
        printf("Process rank: %d; Thread rank: %d\n", rank, omp_get_thread_num());
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
