# Assignment 5, due November 13th 2019

The goal of this assignment is to implement an n-body-type (lagrangian) simulation and consider its optimization potential.

## Exercise 1

This exercise consists in implementing a sequential 2D n-body simulation.

### Description

N-body simulations form a large class of scientific applications, as they are used in research ranging from astrophysics to molecular dynamics. At their core, they model and simulate the interaction of moving particles in physical space. For this assignment, the specific n-body setting relates to astrophysics, where the mutual graviational effect of stars is investigated. Each particle has several properties which include at least
- position,
- velocity, and
- mass.

For each timestep (you can assume `dt = 1`), particles must be moved by first computing the force excerted on them according to the [Newtonian equation for gravity](https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation), `force = G * (mass_1 * mass_2) / radius^2` where `G` is the gravitational constant (and can be assumed as `G = 1` for simplicity). Second, using the computed force on a particle, its position and velocity can be updated via `velocity = velocity + force / mass` and `position = position + velocity`.

### Tasks

- _Provide a sequential implementation of the n-body simulation in 2D. Hints on how to proceed (not mandatory to follow):_
    1. _generate particles randomly, e.g. uniformly distributed_
    2. _provide a function for computing forces and moving particles_
        * Due to Newton's third axiom, the forces between two different particles i and j are given by
          <a href="https://www.codecogs.com/eqnedit.php?latex=\textit{\textbf{F}}_{ij}&space;=&space;-&space;\textit{\textbf{F}}_{ji}&space;=&space;-&space;\frac{m_i&space;\cdot&space;m_j}{r_{ij}^2}&space;\cdot&space;\textit{\textbf{e}}_{\textit{\textbf{r}}_{ij}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\textit{\textbf{F}}_{ij}&space;=&space;-&space;\textit{\textbf{F}}_{ji}&space;=&space;-&space;\frac{m_i&space;\cdot&space;m_j}{r_{ij}^2}&space;\cdot&space;\textit{\textbf{e}}_{\textit{\textbf{r}}_{ij}}" title="\textit{\textbf{F}}_{ij} = - \textit{\textbf{F}}_{ji} = - \frac{m_i \cdot m_j}{r_{ij}^2} \cdot \textit{\textbf{e}}_{\textit{\textbf{r}}_{ij}}" /></a>
        * Together with no particle's mass influencing itself, this means only the forces for i<j or i>j have to be calculated.
        * All vectors have to split into their x & y components. The unit vector from above is therefore given by
          <a href="https://www.codecogs.com/eqnedit.php?latex=\textit{\textbf{e}}_{\textit{\textbf{r}}_{ij}}&space;=&space;\frac{\textit{\textbf{r}}_{ij}}{r_{ij}}&space;=&space;\frac{\left(&space;\left(x_i-x_j\right)&space;\cdot&space;\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;\left(y_i-y_j\right)&space;\cdot&space;\textit{\textbf{e}}_{\text{y}}\right)}{\sqrt{\left(x_i-x_j\right)^2&plus;\left(y_i-y_j\right)^2}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\textit{\textbf{e}}_{\textit{\textbf{r}}_{ij}}&space;=&space;\frac{\textit{\textbf{r}}_{ij}}{r_{ij}}&space;=&space;\frac{\left(&space;\left(x_i-x_j\right)&space;\cdot&space;\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;\left(y_i-y_j\right)&space;\cdot&space;\textit{\textbf{e}}_{\text{y}}\right)}{\sqrt{\left(x_i-x_j\right)^2&plus;\left(y_i-y_j\right)^2}}" title="\textit{\textbf{e}}_{\textit{\textbf{r}}_{ij}} = \frac{\textit{\textbf{r}}_{ij}}{r_{ij}} = \frac{\left( \left(x_i-x_j\right) \cdot \textit{\textbf{e}}_{\text{x}} + \left(y_i-y_j\right) \cdot \textit{\textbf{e}}_{\text{y}}\right)}{\sqrt{\left(x_i-x_j\right)^2+\left(y_i-y_j\right)^2}}" /></a>
        * The force excerted on a particle i is given by                                                 
          <a href="https://www.codecogs.com/eqnedit.php?latex=\textit{\textbf{F}}_{i}&space;=&space;F_{i,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;F_{i,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}&space;=&space;\sum_j{\textit{\textbf{F}}_{ij}}&space;=&space;\sum_j{\left(F_{ij,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;F_{ij,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}\right)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\textit{\textbf{F}}_{i}&space;=&space;F_{i,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;F_{i,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}&space;=&space;\sum_j{\textit{\textbf{F}}_{ij}}&space;=&space;\sum_j{\left(F_{ij,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;F_{ij,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}\right)}" title="\textit{\textbf{F}}_{i} = F_{i, \text{x}}\cdot\textit{\textbf{e}}_{\text{x}} + F_{i, \text{y}}\cdot\textit{\textbf{e}}_{\text{y}} = \sum_j{\textit{\textbf{F}}_{ij}} = \sum_j{\left(F_{ij, \text{x}}\cdot\textit{\textbf{e}}_{\text{x}} + F_{ij, \text{y}}\cdot\textit{\textbf{e}}_{\text{y}}\right)}" /></a>
    3. _move particles in a time loop for a given number of steps_ 
        * The resulting velocity after a time dt=1 is                                      
          <a href="https://www.codecogs.com/eqnedit.php?latex=\textit{\textbf{v}}_{i}&space;=&space;v_{i,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;v_{i,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}&space;=&space;\frac{1}{m_i}\cdot\left(F_{i,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;F_{i,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\textit{\textbf{v}}_{i}&space;=&space;v_{i,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;v_{i,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}&space;=&space;\frac{1}{m_i}\cdot\left(F_{i,&space;\text{x}}\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;F_{i,&space;\text{y}}\cdot\textit{\textbf{e}}_{\text{y}}\right)" title="\textit{\textbf{v}}_{i} = v_{i, \text{x}}\cdot\textit{\textbf{e}}_{\text{x}} + v_{i, \text{y}}\cdot\textit{\textbf{e}}_{\text{y}} = \frac{1}{m_i}\cdot\left(F_{i, \text{x}}\cdot\textit{\textbf{e}}_{\text{x}} + F_{i, \text{y}}\cdot\textit{\textbf{e}}_{\text{y}}\right)" /></a>
        * The resulting position after a time dt=1 is                                   
          <a href="https://www.codecogs.com/eqnedit.php?latex=\textit{\textbf{r}}_{i}&space;=&space;x_i&space;\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;y_i&space;\cdot\textit{\textbf{e}}_{\text{y}}&space;=&space;\left(x_i&plus;v_{i,\text{x}}\right)&space;\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;\left(y_i&plus;v_{i,\text{y}}\right)&space;\cdot\textit{\textbf{e}}_{\text{y}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\textit{\textbf{r}}_{i}&space;=&space;x_i&space;\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;y_i&space;\cdot\textit{\textbf{e}}_{\text{y}}&space;=&space;\left(x_i&plus;v_{i,\text{x}}\right)&space;\cdot\textit{\textbf{e}}_{\text{x}}&space;&plus;&space;\left(y_i&plus;v_{i,\text{y}}\right)&space;\cdot\textit{\textbf{e}}_{\text{y}}" title="\textit{\textbf{r}}_{i} = x_i \cdot\textit{\textbf{e}}_{\text{x}} + y_i \cdot\textit{\textbf{e}}_{\text{y}} = \left(x_i+v_{i,\text{x}}\right) \cdot\textit{\textbf{e}}_{\text{x}} + \left(y_i+v_{i,\text{y}}\right) \cdot\textit{\textbf{e}}_{\text{y}}" /></a>
 
- Measure the execution time for various problem sizes. What can you observe?

## Exercise 2

This exercise consists in investigating and planning optimization and parallelization techniques for Exercise 1.

### Tasks

- Study the nature of the problem in Exercise 1, focusing on its characteristics with regard to optimization and parallelization.
- What optimization methods can you come up with in order to improve the performance of Exercise 1?
- What parallelization strategies would you consider for Exercise 1 and why?

## General Notes

All the material required by the tasks above (e.g. code, figures, etc...) must be part of the solution that is handed in. Your experiments should be reproducible and comparable to your own measurements using the solution materials that you hand in. For source code, please provide a makefile or other, intuitive means of compiling with the required flags and settings.

**Every** member of your group must be able to explain the given problem, your solution, and possible findings. You may also need to answer detailed questions about any of these aspects.

**Please run any benchmarks or heavy CPU loads only on the compute nodes, not on the login node.**
If you want to do some interactive experimentation, use an *interactive job* as outlined in the tutorial. Make sure to stop any interactive jobs once you are done.
