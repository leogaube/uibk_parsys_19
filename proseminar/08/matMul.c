#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define IDX_MATRIX(x,y,nx) ((x)+(y)*(nx))

void initMatrix(int* mat, int ny, int nx);
void transpose(int* mat, int ny, int nx, int* mat_t);
void multiply(int* A, int* B, int*C, int L, int M, int N);
void print_matrix(int* mat, int ny, int nx);


int main(int argc, char **argv) {
    double start = omp_get_wtime();
    int L, M, N = 10;
    if (argc == 2) {
    	L = atoi(argv[1]);
    	M = L;
    	N = L;
    }
    if (argc == 4) {
        L = atoi(argv[1]);
        M = atoi(argv[2]);
        N = atoi(argv[3]);
    }

    int *A = malloc(L*M*sizeof(int));
    int *B = malloc(M*N*sizeof(int));
    int *C = malloc(L*N*sizeof(int)); //result

    srand(1234);
    /*
     *  The shared dimension has to be the fast one.
     *  The typical notation is (#rows,#columns). Therefore, (y,x)
     *  is used for the initialization, keeping in mind that x
     *  (representing a row) is the fast changing index.
     */
    initMatrix(A, L, M);
    initMatrix(B, M, N);

    multiply(A, B, C, L, M, N);

    double end = omp_get_wtime();
    printf("The process took %f seconds to finish. \n", (end - start));

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
	for(int y=0;y<ny;y++){
		for(int x=0;x<nx;x++){
			mat[IDX_MATRIX(x,y,nx)] = (int) ((((double) rand()/RAND_MAX)-0.5)*10);
		}
	}
}


void transpose(int* mat, int ny, int nx, int* mat_t){
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

	// init C
	for(int i=0; i<L*N; i++){
		C[i] = 0;
	}

	// perform the multiplication
	for(int l=0; l<L; l++){
		for(int n=0; n<N; n++){
			for(int m=0; m<M; m++){
				C[IDX_MATRIX(n,l,N)] += A[IDX_MATRIX(m,l,M)]*B_t[IDX_MATRIX(m,n,M)];
			}
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
