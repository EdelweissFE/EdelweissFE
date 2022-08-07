# EdelweissFE: A light-weight, platform-independent, parallel finite element framework.

See the [documentation](https://edelweissfe.github.io/edelweissfe).

EdelweissFE aims at an easy to understand, yet efficient implementation of the finite element method.
Some features are:

 * Python for non performance-critical routines
 * Cython for performance-critical routines
 * Parallelization 
 * Modular system, which is easy to extend
 * Output to Paraview, Ensight, CSV, matplotlib
 * Interfaces to powerful direct and iterative linear solvers

EdelweissFE makes use of the [Marmot](https://github.com/MAteRialMOdelingToolbox/Marmot/) library for finite element and constitutive model formulations.
