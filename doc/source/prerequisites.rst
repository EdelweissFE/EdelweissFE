Prerequisites
#############

EdelweissFE requires 

* Python 3.5+
* NumPy
* SciPy
* SymPy
* Cython 0.28+
* Matplotlib*
* Rich
* Gstools
* OpenMP
* Intel MKL for the PARDISO  (binaries and header files; if you use Anaconda, install packages ``mkl`` and ``mkl-include``)
* `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_
* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ library for linear algebra 

All those packages can be easily installed using Anaconda.
A recent version of Anaconda `<https://anaconda.org/>`_ is sufficient to compile and run EdelweissFE on a Linux system.
Anaconda also provides free access to Intel MKL binaries (via package ``mkl``, usually installed by default) and header files (via package ``mkl-include``, available via ``conda install mkl-include``), so no standalaone installation of the Intel MKL is required.
