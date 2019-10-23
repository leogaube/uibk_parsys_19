# Solution: Exercise 2

## Task 1
*A sequential implementation of a 1-D heat stencil is available in heat_stencil_1D_seq.c. Read the code and make sure you understand what happens. See the Wikipedia article on Stencil Codes for more information.*

Making a calculation with neighboring cells in every time step.

## Task 2
*Consider a parallelization strategy using MPI. Which communication pattern(s) would you choose and why?*

- Splitting the room into sub-rooms requires only little communication between ranks. 
- For a single timestep, most cells can be computed with knowlegde of a single subroom (exception: edge values)
- Each sub-room has to send its edge values to its neighboring sub-rooms using Point-to-Point communication.
- For simplicity, sub-rooms comprising a wall send the edge values closest to the wall to themselves.

*Are there additional changes required in the code beyond calling MPI functions? If so, elaborate!*

- Each rank needs to allocate memory for the subrooms.
- In order to output the temperatures, the entire room needs to be reassembled using `MPI_Gather`.

## Task 3
*Implement your chosen parallelization strategy as a second application heat_stencil_1D_mpi. Run it with varying numbers of ranks and problem sizes and verify its correctness by comparing the output to heat_stencil_1D_seq.*

See programs and outputs. Equality of outputs verified by comparison with sequential program.
If you are insterested in runtime for different values of N=2000-16000, use:
`grep seconds [output_file]`

## Task 4
*Discuss the effects and implications of your parallelization.*

- The `MPI_Gather` to plot the temperature values after 1000 time step takes time. Optimization might be possible if only the final result matters.
- For small problem sizes like N = 2000, parallelization does not make much sense, because there is too much MPI initialization and communication overhead.
- However we have done some experiments for bigger N until N=16000 and yielded a speedup of 2.24 (compared to 1.27 for N=2000) for 4 ranks running on the same cpu.
- Unfortunately the lcc2 cluster was always fully utilized in the last two days; thus we were unable to get reliable results for different number of ranks or even node to node communication!
- Tables & Figures for our results (seq + 4_ranks_same_cpu) can be found in *results.xlsx*
