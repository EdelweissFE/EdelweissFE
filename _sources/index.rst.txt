.. EdelweissFE documentation master file, created by
   sphinx-quickstart on Wed Jul 27 10:10:22 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EdelweissFE
======================

EdelweissFE is a light-weight, platform-independent, parallel finite element framework.
By default, it makes use of the `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_ library for finite element and constitutive model formulations.

EdelweissFE aims to be...

 - ... a development and learning environment for constitutive models and finite elements,
 - ... an easy to use tool for coupled problems,
 - ... a learning platform for learning the finite element method,
 - ... a very flexible tool for implementing and employing special techniques (e.g., the indirect displacement control technique),
   which are often more difficult and time consuming to implement in mature, MPI-parallelized codes,
 - ... an efficient tool for nonlinear simulations up to medium sized problems (~ :math:`10^5` degrees of freedom).

EdelweissFE does not want to ...

 - ... compete with more mature (MPI-parallelized) codes such as `MOOSE <https://mooseframework.inl.gov/>`_,
   which are also compatible with the `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_ library, but less suitable as development environments.
   Usually, EdelweissFE is used for developing and debugging constitutive models until they meet the requirements for production runs on HPC system using more mature frameworks.


.. figure:: borehole.png
   :align: center
   :width: 240px

.. toctree::
   :maxdepth: 2
   :hidden:

   features
   prerequisites
   installation
   documentation/index
   contributors
   publications

Execute a simulation
********************

Run a simulation simply by calling

.. code-block:: console

    python edelweiss.py inputfile.inp

Example input file
******************

.. literalinclude:: ../../examples/CantileverBeam/test3D.inp
   :language: edelweiss
   :caption: File: ``examples/CantileverBeam/test3D.inp``

More example files can be found in the ``testfiles`` folder!

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
