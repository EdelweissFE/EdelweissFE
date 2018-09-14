EdelweissFE: A light-weight, platform-independent, parallel finite element framework.

## Prerequisites

EdelweissFE requires 
- Python 3.5+
- Numpy
- Scipy
- Cython 0.28+
- Intel MKL for the PARDISO solver
- Matplotlib
- OpenMP 
- bftUserLibrary (Elements and Materials)

A recent version of Anaconda (https://anaconda.org/) is sufficient to compile and run EdelweissFE on a linux system.

## Configuration

Modify install.py:

- define the parent directory containing the `bftUserLibrary` directory (e.g., `/home/user/constitutiveModelling`)

## Installation

EdelweissFE depends on several Cython modules, which must be compiled prior to a first run.

Simply run
`python install.py`


## Run the validation examples

`python validateEdelweiss.py`

## Execute a simulation

`python edelweiss.py INPUT.inp`

# Parallelization

If a parallel solver (e.g, NISTParallel, NISTPArcLength) is selected in the .inp file, EdelweissFE  automatically determines the max. number of threads dependent on the host architecture.
However, it is RECOMMENDED to enforce a fixed number of threads by running

`OMP_NUM_THREADS=XX python edelweiss.py INPUT.inp`

This ensures that the same number of threads XX is employed both in EdelweissFE as well as in the underlying Intel MKL (e.g., if the PARDISO linear solver is used).

