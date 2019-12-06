#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

// Function to check if two queens threaten each other or not
int is_safe(int N, int queen_rows[N]) {
	for (int i = 0; i < N; i++)
	{
		for (int j = i+1; j < N; j++)
		{
			// two queens in the same row => not a solution!
			if (queen_rows[i] == queen_rows[j]) return 0;
			
			// two queens in the same diagonal => not a solution!
			if (queen_rows[i] - queen_rows[j] == i - j ||
			    queen_rows[i] - queen_rows[j] == j - i)
			    return 0;
		}
	}

	return 1;
}

int main(int argc, char* argv[])
{
    int N = 8;
    if (argc == 2) {
      N = atoi(argv[1]);
    }
    int max_iter = 1;
    
    clock_t start = clock();
    int solutions = 0;
	        
    for (int i = 0; i < N; i++)
    {
        max_iter *= N;
    }

	for (int iter = 0; iter < max_iter; iter++)
	{
		int code = iter;
		int i;
	    int queen_rows[N];
		// the index correspond to the queen's number and the queen's collumn
		// we only generate configurations where there's only one queen per collumn
		for (i = 0; i < N; i++)
		{
			queen_rows[i] = code % N;
			
			code /= N;
		}
		
		if (is_safe(N, queen_rows))
		{
            solutions++;
		}
	}

    clock_t end = clock();
    printf("The result is: %d\n", solutions);
    printf("The process took %f seconds to finish. \n", ((double)(end - start)) / CLOCKS_PER_SEC);    
    
	return 0;
}
