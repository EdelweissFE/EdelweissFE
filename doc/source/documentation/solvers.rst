Solvers
===========

.. automodule:: fe.config.solvers
   :members: 

.. pprint:: fe.config.solvers.solverLibrary

Choose the solver in the ``*job`` definition:

.. code-block:: console

    *job, name=myJob, domain=2d, solver=NISTParallel


``NIST`` - Nonlinear Implicit Static
------------------------------------

.. automodule:: fe.solvers.nonlinearimplicitstatic
   :members: 


``NISTParallelForMarmotElements`` - Nonlinear Implicit Static (parallel)
------------------------------------------------------------------------

.. automodule:: fe.solvers.nonlinearimplicitstaticparallelmk2
   :members: 

``NISTPArcLength`` - Nonlinear Implicit Static - Arc length 
-----------------------------------------------------------

.. automodule:: fe.solvers.nonlinearimplicitstaticparallelarclength
   :members: 
