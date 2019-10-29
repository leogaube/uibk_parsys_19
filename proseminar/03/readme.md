# Assignment 3, due October 30th 2019

The goal of this assignment is to extend the heat stencil application and measure its performance.

## Exercise 1

This exercise consists in extending the heat stencil application of Assignment 2 to two and three dimensions.

### Tasks

- Extend the heat stencil application to the two- and three-dimensional cases and name them `heat_stencil_2D` and `heat_stencil_3D`.
- Provide sequential and MPI implementations and run them with multiple problem and machine sizes.
- How can you verify the correctness of your applications?

## Solution 1

## Task 1&2 (Implementation)
- See mpi_2D/3D and seq_2D/3D.

## Task 3 (Verification)
Idea: Use already verified 1D version to verify in more dimensions.

### Verification along one axis
Running the 2D code with the dimensions (Nx, 1) and the 3D code with (Nx, 1, 1) should result
in the same output as the 1D code.

Results:
* see `out_verification2D.dat` and `out_verification3D.dat`
* The results match well as long as the heat is not **hitting a wall**.
    * The assumption that the wall has the same temperature as the first/last cell of the room 
leads to a reflection of the heat towards the center and therefore a deviation.
    * If this is not desired/valid, **other boundary conditions** would have to be chosen.
    * Periodic boundary conditions help, but lead to multiple heat sources due to periodicity.

### Verification in multiple dimensions
Given the program is valid, a point source should produce a spherically symmetric heat distribution. 
Therefore, calculating the problem in 1D and 'rotating' the result with the origin as anchor should result
in the same distribution.

Results:
- The same problems as in the verification along one axis.
- Furthermore it seems like the function T = Tc + 0.2 ( Tl + Tr + Tu + Td + Tf + Tb + (-6 Tc )) is **not valid**:
    - If the difference between Tc and the surrounding cells is too large, the central cell could cool down
    below the T of the coldest neighboring cell. If all neighboring cells are at T=0, the function yields
    T = - 0.2 Tc. -> **Negative temperatures** are not physical.
    - The value 0.2 or the whole function have to be adapted.
- In 2D this effect is not existing as T = Tc + 0.2 ( Tl + Tr + Tu + Td + (-4 Tc )) can at most evaluate
to T=0.2Tc. An overestimation of the cooling/heating effect of the neighboring cells might however still occur.


## Exercise 2

This exercise consists in measuring all heat stencil variants (1D, 2D and 3D) to get a grasp of their performance behavior.

### Tasks

- Measure the speedup and efficiency of all three stencil codes for varying problem and machine sizes/mappings. Consider using strong scalability, weak scalability, or both. Justify your choice.
- Illustrate the data in appropriate figures and discuss them. What can you observe?
- Bonus question: Measure and illustrate an application throughput metric. What can you observe?

## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an *interactive job* as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
