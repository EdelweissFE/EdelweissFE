# EdelweissFE: A light-weight, platform-independent, parallel finite element framework.

EdelweissFE makes use of the [Marmot](https://github.com/MAteRialMOdelingToolbox/Marmot/) library for finite element and constitutive model formulations.

## Prerequisites

EdelweissFE requires 
- Python 3.5+*
- NumPy*
- SciPy*
- SymPy*
- Cython 0.28+*
- Matplotlib*
- Rich
- OpenMP*
- Intel MKL for the PARDISO*  (binaries and header files; if you use Anaconda, install packages `mkl` and `mkl-include`)
- [Marmot](https://github.com/MAteRialMOdelingToolbox/Marmot/)
- Eigen Library for Linear Algebra (http://eigen.tuxfamily.org/index.php?title=Main_Page)

*: Provided via Anaconda packages
A recent version of Anaconda (https://anaconda.org/) is sufficient to compile and run EdelweissFE on a Linux system.
Anaconda also provides free access to Intel MKL binaries (via package `mkl`, usually installed by default) and header files (via package `mkl-include`, available via ```conda install mkl-include```), so no standalaone installation of the Intel MKL is required.

## Configuration

Customize setup.py by defining all paths pointing to the respective libraries.
Default paths are already defined, and usually only minor modifications should be required.

## Installation

EdelweissFE depends on several Cython modules, which must be compiled prior to first running the software.

Simply run
`python setup.py build_ext -i`

Enforce recomipilation and installation with
`python setup.py build_ext -i --force`

## Run the validation examples

`python validateEdelweiss.py`

Recreate the validation reference solutions (please only if you know what you are doing)

`python validateEdelweiss.py --create`

## Execute a simulation

`python edelweiss.py INPUT.inp`

### Parallelization

If a parallel solver (e.g, NISTParallel, NISTPArcLength) is selected in the .inp file, EdelweissFE  automatically determines the max. number of threads dependent on the host architecture.
However, it is RECOMMENDED to enforce a fixed number of threads by running

`OMP_NUM_THREADS=XX python edelweiss.py INPUT.inp`

This ensures that the same number of threads XX is employed both in EdelweissFE as well as in the underlying Intel MKL (e.g., if the PARDISO linear solver is used).

