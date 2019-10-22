# Solution: Assignment 2

## Task 1
A sequential implementation of a 1-D heat stencil is available in heat_stencil_1D_seq.c. Read the code and make sure you understand what happens. See the Wikipedia article on Stencil Codes for more information.

Making a calculation with neighboring cells in every time step.

## Task 2
Consider a parallelization strategy using MPI. Which communication pattern(s) would you choose and why? Are there additional changes required in the code beyond calling MPI functions? If so, elaborate!

- Splitting the room into sub-rooms.
- Sub-rooms have to send the edge values to the neighboring sub-rooms.
- Sub-rooms on the edges send themselves the edge values for simplicity.
- After every time step the complete room is build up to do the plotting.

## Task 3
Implement your chosen parallelization strategy as a second application heat_stencil_1D_mpi. Run it with varying numbers of ranks and problem sizes and verify its correctness by comparing the output to heat_stencil_1D_seq.
See programs and outputs. Equality of outputs verified by comparison with separate program.

## Task 4
Discuss the effects and implications of your parallelization.

- The room size has to be divisible by the number of ranks.
- The gather to plot the values after each time step takes time. Optimization possible if not every time step has to be plotted. 
