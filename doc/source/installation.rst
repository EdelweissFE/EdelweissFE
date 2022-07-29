Installation
============

Marmot
******
By default, a working instance of `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_, 
containing specific implementations of finite elements and constutitive models,
is required.

Please build Marmot before you build EdelweissFE.

Configuration
*************

Customize setup.py by defining all paths pointing to the respective libraries.
Default paths are already defined, and usually only minor modifications should be required.

Building EdelweissFE
********************

EdelweissFE depends on several Cython modules, which must be compiled prior to first running the software.

Simply run

.. code-block:: console

    python setup.py build_ext -i

Enforce recompilation and installation with

.. code-block:: console

    python setup.py build_ext -i --force

Run the validation examples

.. code-block:: console

    python run_tests.py

Recreate the validation reference solutions (please only if you know what you are doing)

.. code-block:: console

    python run_tests.py --create


TLDR
****

Assuming that you are in an empty directory,
you can quickly get a working version of EdelweissFE in a Linux based 
environment:

Steps
_____

Get mamba:

.. code-block:: console

    curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
    bash Mambaforge-Linux-x86_64.sh -b -p ./mambaforge3

Add mamba to your environment (repeat this step if you close your shell):

.. code-block:: console

    export PATH=$PWD/mambaforge3/bin:$PATH

Install necessary mamba packages:

.. code-block:: console

    mamba install --yes python=3.10 cython numpy sympy matplotlib scipy 
    mamba install --yes compilers cmake gfortran gstools  
    mamba install --yes -c anaconda mkl mkl-include intel-openmp
    mamba install --yes sphinx sphinx_rtd_theme

Get Eigen (for EdelweissFE and Marmot):

.. code-block:: console

    git clone   https://gitlab.com/libeigen/eigen.git
    cd eigen
    mkdir build
    cd build
    cmake -DBUILD_TESTING=OFF  -DINCLUDE_INSTALL_DIR=$(python -c "import sys; print(sys.prefix)")/include -DCMAKE_INSTALL_PREFIX=$(python -c "import sys; print(sys.prefix)") ..
    make install
    cd ../../

Get autodiff (for Marmot):

.. code-block:: console

    git clone  https://github.com/autodiff/autodiff.git
    cd autodiff
    mkdir build
    cd build
    cmake -DAUTODIFF_BUILD_TESTS=OFF -DAUTODIFF_BUILD_PYTHON=OFF -DAUTODIFF_BUILD_EXAMPLES=OFF -DAUTODIFF_BUILD_DOCS=OFF -DCMAKE_INSTALL_PREFIX=$(python -c "import sys; print(sys.prefix)") ..
    make install
    cd ../

.. Get Fastor:

.. .. code-block:: console

..     git clone https://github.com/romeric/Fastor.git
..     cd Fastor
..     cmake -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=$(python -c "import sys; print(sys.prefix)") .
..     make install
..     cd ../

Get Marmot: 

.. code-block:: console

    git clone --recurse https://github.com/MAteRialMOdelingToolbox/Marmot.git
    cd Marmot
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=$(python -c "import sys; print(sys.prefix)") ..
    make install
    cd ../../

Get, build and test EdelweissFE:

.. code-block:: console

    git clone https://github.com/EdelweissFE/EdelweissFE.git
    cd EdelweissFE
    python setup.py build_ext -i
    python run_tests.py

Build this documentation:

.. code-block:: console

    sphinx-build ./doc/source/ ./docs -b html
