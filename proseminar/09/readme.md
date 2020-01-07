# Assignment 9, due December 11th 2019

The goal of this assignment is to implement the n-queens problem using OpenMP.

## Exercise 1

### Description

N-queens is a popular branch-and-bound problem. The goal is to compute the number of possible ways to place N chess queens on a NxN chess board without them attacking each other. See https://en.wikipedia.org/wiki/Eight_queens_puzzle for further information.

### Tasks

- **Implement a sequential version of the n-queens problem. Benchmark your program for several problem sizes. What can you observe?**

Runtime increases exponensionally with the problem size.
A 16-queens problem takes about 9 minutes, whereas a 10-queens problem only takes about 10 milliseconds.

- **Parallelize your program using OpenMP. Which OpenMP constructs are suitable candidates and why?**

The "openmp for" and "openmp tasks" directive are both suitable candidates, but we chose tasks, because they support easy nested  work sharing (e.g. for a backtracking problem like n-queens)!

For each task, we had to use copies of the chess board in order to use them truly in parallel!
Unfortunately `firstprivate(board)` does not seem to be enough, because the thread-private board variable is only initialized once the task starts (?), so the global board variable might have changed by that time.

We were unable to figure out how to bypass this problem, so we solve n-queens in parallel by creating N different boards in the beginning each with a queen in a different position in the first row.

- **Benchmark your problem for several problem sizes and numbers of threads. What can you observe?**

Due to our specific setup, our parallel program has good efficiency if N is divisable by the number of ranks. That way there is low load imbalance and we can reach an efficiency of 98% on a 16x16 board with 8 ranks.
In contrast efficiency for the same problem size with only 6 ranks is worse (90%) due to more load imbalance.

We assume this could be solved with more Tasks (e.g. not just N, but N^2) by simply declaring a task with `final(row >= 2)` for each resursive call of n_queens, but - as mentioned above - our efforts failed at dynamically initiallizing thread-private copies of board with `firstprivate(board)`!

We also get really consistent runtime results on all 6 runs, which reaffirm out assumptions. 

Our plots can be found in the `./results` directory.

- **What are potential optimizations for this application?**

We have already done some optimizations in our sequential program, by only adding a single queen at a time for each row and also by checking whether a new queen would be valid for a subboard NxR (with R=row of new queen), instead of the entire NxN board.

I could also be possible to use symmetry to solve n-queens faster, but that might be hard to implement using backtracking!

The german Wikipedia articel also mentions an _iterative repair algorithm_ (starting from an invalid position). However, that algorithm only makes sense for single (non-trivial) solutions of large n-queen problems and it is not even garanteed to always find a solution!

## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an *interactive job* as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
