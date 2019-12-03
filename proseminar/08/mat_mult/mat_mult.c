#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

#define IDX_MATRIX(x,y,nx) ((x)+(y)*(nx))

void initMatrix(int* mat, int ny, int nx);
void transpose(int* mat, int ny, int nx, int* mat_t);
void multiply(int* A, int* B, int*C, int L, int M, int N);
void print_matrix(int* mat, int ny, int nx);


int main(int argc, char **argv) {
#ifdef OMP
    double start = omp_get_wtime();
#else
    clock_t start = clock();
#endif
    int L, M, N;
    L = M = N = 10;
    if (argc >= 2) {
    	L = M = N = atoi(argv[1]);
    }
    if (argc >= 4) {
        L = atoi(argv[1]);
        M = atoi(argv[2]);
        N = atoi(argv[3]);
    }
    printf("L=%d, M=%d, N=%d\n", L,M,N);

    int *A = malloc(L*M*sizeof(int));
    int *B = malloc(M*N*sizeof(int));
    int *C = malloc(L*N*sizeof(int)); //result

	/*
	 *  The shared dimension has to be the fast one.
	 *  The typical notation is (#rows,#columns). Therefore, (y,x)
	 *  is used for the initialization, keeping in mind that x
	 *  (representing a row) is the fast changing index.
	 */
	initMatrix(A, L, M);
	initMatrix(B, M, N);

	multiply(A, B, C, L, M, N);

#ifdef OMP
    double end = omp_get_wtime();
    printf("The process took %f seconds to finish. \n", (end - start));
#else
    clock_t end = clock();
    printf("The process took %f seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);
#endif

#ifdef VERBOSE
    printf("A:\n");
    print_matrix(A, L, M);
    printf("times B:\n");
    print_matrix(B, M, N);
    printf("gives C:\n");
    print_matrix(C, L, N);
#endif
    free(A);
    free(B);
    free(C);
    return EXIT_SUCCESS;
}


void initMatrix(int* mat, int ny, int nx){
#ifdef OMP
	unsigned int seeds[8] = {1, 2, 3, 4, 5, 6, 7, 8}; //up to 8 threads are possible
#else
	srand(1234);
#endif
	#pragma omp parallel for
	for(int y=0;y<ny;y++){
		for(int x=0;x<nx;x++){
#ifdef OMP
			double rnd = (double)rand_r(&seeds[omp_get_thread_num()]);
#else
			double rnd = (double)rand();
#endif
			mat[IDX_MATRIX(x,y,nx)] = (int) (((rnd/RAND_MAX)-0.5)*10);
		}
	}
}


void transpose(int* mat, int ny, int nx, int* mat_t){
	#pragma omp parallel for
	for(int y=0; y<ny; y++){
		for(int x=0; x<nx; x++){
			mat_t[IDX_MATRIX(y,x,ny)] = mat[IDX_MATRIX(x,y,nx)];
		}
	}
}

/**
 * assuming the dimensions
 * A = (L,M)
 * B = (M,N)
 * C = (L,N)
 */
void multiply(int* A, int* B, int* C, int L, int M, int N){
	//swap columns and rows of B for the cash as the second index (x) is the fast changing one
	int* B_t = malloc(N*M*sizeof(int));
	transpose(B, M, N, B_t);

	// perform the multiplication
	#pragma omp parallel for
	for(int l=0; l<L; l++){
		for(int n=0; n<N; n++){
			int sum=0;
			for(int m=0; m<M; m++){
				sum += A[IDX_MATRIX(m,l,M)]*B_t[IDX_MATRIX(m,n,M)];
			}
			C[IDX_MATRIX(n,l,N)] = sum;
		}
	}

	free(B_t);
}


void print_matrix(int* mat, int ny, int nx){
	printf("-----------------------------------------\n");
	for(int y=0; y<ny; y++){
		for(int x=0; x<nx; x++){
			printf("%d\t",mat[IDX_MATRIX(x,y,nx)]);
		}
		printf("\n");
	}
	printf("-----------------------------------------\n");
}
