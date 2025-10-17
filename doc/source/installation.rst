Installation
============

Marmot
******
By default, a working instance of `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_,
containing specific implementations of finite elements and constitutive models,
is required.
Marmot itself requires the `Eigen <https://eigen.tuxfamily.org/>`_ library,
and potentially `Fastor <https://github.com/romeric/Fastor>`_, depending on the requested modules.

Please build Marmot before you build EdelweissFE.

Configuration
*************

Customize ``setup.py`` by defining all paths pointing to the respective libraries.
Default paths are already defined, and usually only minor modifications should be required.

Building EdelweissFE
********************

EdelweissFE depends on several Cython modules, which must be compiled prior to running the EdelweissFE.

To build EdelweissFE and install it using pip, simply run from within the main folder:

.. code-block:: console

    cd ./EdelweissFE
    pip install .

Run a simulation with an inputfile using

.. code-block:: console

    edelweissfe your_input_file.inp

Run all the validation examples

.. code-block:: console

    run_tests_edelweissfe ./testfiles/

Recreate the validation reference solutions (only if you know what you are doing)

.. code-block:: console

    run_tests_edelweissfe --create

Alternatively, to build the modules (inplace) and not install them using pip, simply run

.. code-block:: console

    python setup.py build_ext -i

Enforce a recompilation with

.. code-block:: console

    python setup.py build_ext -i --force


TLDR
****

Assuming that you are in an empty directory,
you can quickly get a working version of EdelweissFE in a Linux based
environment:

Installation steps
__________________

If necessary, get `Anaconda <https://www.anaconda.com/>`_

.. code-block:: console
   :caption: Step 1

    curl -L -O \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh
    bash Miniforge3-Linux-aarch64.sh -b -p ./miniforge3

Add mamba to your environment:

.. code-block:: console
   :caption: Step 2

    export EWROOT=$PWD
    export PATH=$EWROOT/mambaforge3/bin:$PATH
    conda init --all
    exit

Restart shell and activate mamba

.. code-block:: console
   :caption: Step 3

    export EWROOT=$PWD
    conda activate

Get EdelweissFE:

.. code-block:: console
   :caption: Step 4

    git clone https://github.com/EdelweissFE/EdelweissFE.git

Install necessary mamba packages:

.. code-block:: console
   :caption: Step 5

    mamba install --file EdelweissFE/requirements.txt

Get Eigen (for EdelweissFE and Marmot):

.. code-block:: console
   :caption: Step 6

    cd $EWROOT
    git clone --branch 3.4.0  https://gitlab.com/libeigen/eigen.git
    cd eigen
    mkdir build
    cd build
    cmake \
        -DBUILD_TESTING=OFF  \
        -DINCLUDE_INSTALL_DIR=$CONDA_PREFIX/include \
        -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
        ..
    make install

Get autodiff (for Marmot):

.. code-block:: console
   :caption: Step 7

    cd $EWROOT
    git clone --branch v1.1.0 https://github.com/autodiff/autodiff.git
    cd autodiff
    mkdir build
    cd build
    cmake \
        -DAUTODIFF_BUILD_TESTS=OFF \
        -DAUTODIFF_BUILD_PYTHON=OFF \
        -DAUTODIFF_BUILD_EXAMPLES=OFF \
        -DAUTODIFF_BUILD_DOCS=OFF \
        -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
        ..
    make install

Get Fastor:

.. code-block:: console
   :caption: Step 8

    cd $EWROOT
    git clone https://github.com/romeric/Fastor.git
    cd Fastor
    cmake -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX .
    make install
    cd ../

Get Marmot:

.. code-block:: console
   :caption: Step 9

    cd $EWROOT
    git clone --recurse https://github.com/MAteRialMOdelingToolbox/Marmot.git
    cd Marmot
    mkdir build
    cd build
    cmake \
        -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
        ..
    make install

Build and test EdelweissFE:

.. code-block:: console
   :caption: Step 10

    cd $EWROOT
    cd EdelweissFE
    pip install .
    run_tests_edelweissfe ./testfiles/

Build this documentation:

.. code-block:: console
   :caption: Step 11

    sphinx-build ./doc/source/ ./docs -b html
