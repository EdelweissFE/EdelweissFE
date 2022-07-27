.. EdelweissFE documentation master file, created by
   sphinx-quickstart on Wed Jul 27 10:10:22 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EdelweissFE !
========================

EdelweissFE is a light-weight, platform-independent, parallel finite element framework.
EdelweissFE makes use of the `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_ library for finite element and constitutive model formulations.

.. toctree::
   :maxdepth: 2
   :hidden:

   prerequisites
   installation
   features
   documentation
   syntax
   contributors

Execute a simulation
********************

Run a simulation simply by calling

.. code-block:: console

    python edelweiss.py INPUT.inp

Example input file
******************


.. literalinclude:: ../../testfiles/LinearElasticIsotropic/test.inp
   :language: console

More example files can be found in the ``testfiles`` folder!

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
