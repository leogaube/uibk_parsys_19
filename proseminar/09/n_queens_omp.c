#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>

_Atomic int global = 0;

// Function to check if two queens threaten each other or not
int is_safe(int N, char board[N][N], int row, int col) {
    // return false if two queens share the same column
    for (int i = 0; i < row; i++) {
        if (board[i][col] == 'Q') {
            return 0;
        }
    }

    // return false if two queens share the same \ diagonal
    for (int i = row, j = col; i >= 0 && j >= 0; i--, j--) {
        if (board[i][j] == 'Q') {
            return 0;
        }
    }

    // return false if two queens share the same / diagonal
    for (int i = row, j = col; i >= 0 && j < N; i--, j++) {
        if (board[i][j] == 'Q') {
            return 0;
        }
    }

    return 1;
}

void set_queen(int N, char board[N][N], int row, int* solutions) {
    // if N queens are placed successfully, add one to solutions
    if (row == N) {
#pragma omp atomic
        *solutions += 1;
        global++;
        return;
    }

    // place Queen at every square in current row r
    // and recur for each valid movement
    for (int i = 0; i < N; i++) {
        // if no two queens threaten each other
        if (is_safe(N, board, row, i)) {
            // place queen on current square
            board[row][i] = 'Q';
            set_queen(N, board, row + 1, solutions);
            board[row][i] = '-';
        }
    }
}

void n_queen(int N, int* solutions) {
    // place Queen at every square in starting row
#pragma omp taskloop shared(solutions)
    for (int i = 0; i < N; i++) {
        char board[N][N];
        memset(board, '-', sizeof board);
        board[0][i] = 'Q';
        set_queen(N, board, 1, solutions);
    }
}

int main(int argc, char** argv) {
    double start = omp_get_wtime();

    // 'parsing' optional input parameter = problem size
    int N = 8;
    if (argc == 2) {
        N = atoi(argv[1]);
    }
    if (argc == 3) {
        N = atoi(argv[1]);
        omp_set_num_threads(atoi(argv[2]));
    }

    int* solutions = malloc(sizeof(int));
    *solutions = 0;

#pragma omp parallel
#pragma omp single
    n_queen(N, solutions);

    double end = omp_get_wtime();
    printf("The result is: %d\n", *solutions);
    printf("The global result is: %d\n", global);
    printf("The process took %f seconds to finish. \n", end - start);

    free(solutions);

    return EXIT_SUCCESS;
}