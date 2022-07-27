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

Enforce recomipilation and installation with

.. code-block:: console

    python setup.py build_ext -i --force

Run the validation examples

.. code-block:: console

    python run_tests.py

Recreate the validation reference solutions (please only if you know what you are doing)

.. code-block:: console

    python run_tests.py --create
