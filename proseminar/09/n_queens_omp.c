#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>

_Atomic int global = 0;

void print_board(int N, char board[N][N]){
    printf("\n");
    for (int row = 0; row < N; row++)
    {
        for (int col = 0; col < N; col++)
        {
            printf("%c", board[row][col]);
        }
        printf("\n");
    }
    printf("\n");
}

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

void n_queens(int N, char board[N][N], int row, int* solutions) {
    // if N queens are placed successfully, add one to solutions
    if (row == N) {
#ifdef VERBOSE
        #pragma omp critical
        print_board(N, board);
#endif

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
            n_queens(N, board, row + 1, solutions);
            board[row][i] = '-';
        }
    }
}

void setup_n_queens(int N, int* solutions) {
    // create N sets for boards each with a single queen placed in the first row and run the n-queens algorithm as tasks
    for (int i = 0; i < N; i++) {
        char board[N][N];
        memset(board, '-', sizeof board);
        board[0][i] = 'Q';
        #pragma omp task
        n_queens(N, board, 1, solutions);
    }
}

int main(int argc, char** argv) {
    double start = omp_get_wtime();

    // 'parsing' optional input parameter = problem size
    int N = 8;
    if (argc == 2) {
        N = atoi(argv[1]);
    }

    int* solutions = malloc(sizeof(int));
    *solutions = 0;

    #pragma omp parallel
    #pragma omp single
    setup_n_queens(N, solutions);

    double end = omp_get_wtime();
    printf("The result is: %d\n", *solutions);
    printf("The global result is: %d\n", global);
    printf("The process took %f seconds to finish. \n", end - start);

    free(solutions);

    return EXIT_SUCCESS;
}