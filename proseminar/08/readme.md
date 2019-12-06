# Assignment 8, due December 4th 2019

The goal of this assignment is to parallelize several applications with OpenMP.

## Exercise 1

### Tasks

- _Use OpenMP to parallelize the Monte Carlo computation of π of Assignment 2 and the 2D heat stencil simulation of Assignment 3._

See the implementations in the monte_carlo and heat_stencil folders.

- _Measure the execution time of your OpenMP programs for several problem sizes and for 1 to 8 threads._

See results.

- _Illustrate the data in appropriate speedup/efficiency figures and discuss them. What can you observe?_

The results are quite inconclusive as for whatever reason, they fluctuate a lot and often give efficiency values above 1. This is very unlikely as it happens even for the most simplistic versions and not in a very reproducible way. We suspect that there is some other optimization in the background that boosts the values, or an error that we have not found/corrected yet.

- _Try to maximize their performance by considering all sequential and parallelism-related optimizations we discussed so far. Which did you choose and why?_

As the data is so inconclusive, a scientific optimization with benchmarking steps in between is not feasible. We performed different optimization steps in the second task however to show how this would be done in principle.

## Exercise 2

### Tasks

- _Use OpenMP to develop a parallel matrix multiplication program._

See the mat_mult folder.

- _Measure the execution time of your OpenMP programs for several matrix sizes and for 1 to 8 threads._

See the results folder.

- _Illustrate the data in appropriate speedup/efficiency figures and discuss them. What can you observe?_

The same inconclusive behavior as above. As a remark, this happens only lcc2. We also have results & figures run on a local machine (up to 4 ranks), which show less flactuation between multiple runs (see `local_results` branch). However, we still only get an efficiency of 84% on such a trivial parallelization, which might be caused by different turbo boost frequencies for varying number of threads (next time we will disable it)!

- _Try to maximize the performance by considering all sequential and parallelism-related optimizations we discussed so far. Which did you choose and why?_

Same problem as before, however several steps of optimization have been performed:

- V1 is the most simple version parallelizing only the outermost for loops.
- V2 includes the collapse instructions.
- V3 does not re-create the parallel environment.
- V4 includes a static scheduling as the workload for each loop-run should be similar.

## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an *interactive job* as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
