#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>

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

void n_queen(int N, char board[N][N], int row, int* solutions) {
	// if N queens are placed successfully, add one to solutions
	if (row == N) {
    *solutions += 1;
		return;
	}

	// place Queen at every square in current row r
	// and recur for each valid movement	
  for (int i = 0; i < N; i++) {
    // if no two queens threaten each other
    if (is_safe(N, board, row, i)) {
      // place queen on current square
      board[row][i] = 'Q';

      n_queen(N, board, row + 1, solutions);

      // backtrack and remove queen from current square
      board[row][i] = '-';
    }
  }
}

int main(int argc, char** argv) {
  clock_t start = clock();

    // 'parsing' optional input parameter = problem size
  int N = 8;
  if (argc == 2) {
    N = atoi(argv[1]);
  }

	// mat[][] keeps track of position of Queens in current configuration
	char board[N][N];
  int *solutions = malloc(sizeof(int));
  *solutions = 0;

	// initalize mat[][] by '-'
	memset(board, '-', sizeof board);

  n_queen(N, board, 0, solutions);


  clock_t end = clock();
  printf("The result is: %d\n", *solutions);
  printf("The process took %f seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);

  free(solutions);

  return EXIT_SUCCESS;
}