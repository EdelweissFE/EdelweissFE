Prerequisites
#############

It is recommended to use a recent version of `Anaconda <https://anaconda.org/>`_ to compile and run EdelweissFE on a Linux system.
Then EdelweissFE requires the packages

.. include:: ../../requirements.txt
   :literal:

Anaconda provides free access to Intel MKL binaries (via package ``mkl``, usually installed by default) and header files (via package ``mkl-include``, available via ``conda install mkl-include``), so no standalaone installation of the Intel MKL is required.

In addition

* `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_
* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ library for linear algebra

are usually required.
