# Assignment 4, due November 6th 2019

The goal of this assignment is to extend the 3D heat stencil application of Assignment 3 to include multiple domain decomposition schemes and investigate the effect of non-blocking communication.

## Exercise 1

### Tasks

- *Provide implementations for the 3D heat stencil application that rely on the three domain decomposition variants presented in the lecture.*

Our source code can be found in `slabs.c` and `cubes.c` each with a dependency to `heat_stencil`. Note that `cubes.c` also allows non-cubic number of ranks and it will automatically determine an appropriate layout (e.g. 16 ranks --> 2\*2\*4). `cubes.c` will also act as source code for `poles`, depending on the exebutable name (poles will set ranks in x-dimension Px=1). 

- *Measure their speedup and efficiency for multiple problem and machine sizes as in the previous exercise.*

Figures can be found in `/proseminar/04/results` as interactive html-graphs. By clicking on a legend-group you can disable/enable it in the graph. For clarity, make sure to only select one or two groups at a time! \
Each figure comprises the same data, but groups them differently (by #ranks, by domain decomposition and for single evaluation).

- *Illustrate the data in appropriate figures and discuss them. What can you observe?*

Our *slabs* algorithm outperforms *cubes/poles* in almost every test (with the exception of really small room sizes and many #ranks e.g. 64x64x64 and 64 ranks). The only explanation we could come up with was that cubes and poles have really bad cache efficiency due to accessing x- and y-slices of each subroom in every iteration for sending data to neighbouring ranks. *slabs* would probably also suffer from cache efficiency if we where to implement it by splitting the room in the x-dimension rather then the z-dimension!

*slabs* even achieves efficiencies of over 95% for room sizes of 256x256x256 and #ranks <= 8, however efficiency rapidly decreases for more ranks (low strong scalability), especially with small room_sizes. 

*cubes* and *poles* on the other hand only achieves an efficiency of about 47% for the same room_size/#ranks, however it strongly scales a lot better (there is only a mild decrease in efficiency for more ranks)

Weak scalability is also quite good for the *slabs* and *poles/cubes* algorithms: bigger room sizes also mean better efficiency. Note that our figures suggest that this trend might not continue for room sizes bigger than 512x512x512 and require more tests!

## Exercise 2

### Tasks

- *Switch from blocking communication to non-blocking communication (or vice versa) and measure the application again.*

Our source code can be found in `slabs_blocking.c` and `cubes_blocking.c` each with a dependency to `heat_stencil`. 

- *Illustrate the data in appropriate figures and discuss them. Can you observe any difference?*

As we would expect, there does not seem to be any significent difference between blocking and non-blocking variants, because we have perfectly balanced loads and we are using the same hardware for each rank. This would probably change in favor of non-blocking communication, if we were to allow room sizes that would result in unbalanced loads!

## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an *interactive job* as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
