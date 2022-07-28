Stepactions
===========

Stepactions are defined within a ``*step`` defintion, 
by specifying their ``name`` and a list of ``option=value``, for example

.. code-block:: console

    *step, jobName=myJob, maxInc=1e-0, minInc=1e-7, maxNumInc=200, maxIter=12, stepLength=1

    dirichlet, name=bottom, nSet=gen_bottom,  field=displacement, 2=0, 1=0
    dirichlet, name=top, nSet=gen_rightTop, field=displacement, 2=0, 

    distributedload, name=dloadt, surface=gen_top, type=pressure, magnitude=0, f(t)=t
    distributedload, name=dloadl, surface=gen_left, type=pressure, magnitude=0, f(t)=t


``bodyforce`` - Bodyforce loads 
-------------------------------

.. automodule:: fe.stepactions.bodyforce
   :members: __doc__

.. pprint:: fe.stepactions.bodyforce.documentation


``dirichlet`` - Standard birichlet BC
-------------------------------------

.. automodule:: fe.stepactions.dirichlet
   :members: __doc__

.. pprint:: fe.stepactions.dirichlet.documentation


``distributedload`` - surfaces loads
------------------------------------

.. automodule:: fe.stepactions.distributedload
    :members: __doc__

.. pprint:: fe.stepactions.distributedload.documentation

``geostatic`` - Geostatic stress definition
-------------------------------------------

.. automodule:: fe.stepactions.geostatic
    :members: __doc__

.. pprint:: fe.stepactions.geostatic.documentation

``indirectcontrol`` - Indirect displacement control
---------------------------------------------------

.. automodule:: fe.stepactions.indirectcontrol
    :members: __doc__

.. pprint:: fe.stepactions.indirectcontrol.documentation

``indirectcontractioncontrol`` - Indirect displacement -- contraction ring control
----------------------------------------------------------------------------------

.. automodule:: fe.stepactions.indirectcontractioncontrol
    :members: __doc__

.. pprint:: fe.stepactions.indirectcontractioncontrol.documentation

``initalizemateiral`` - Initialize materials
--------------------------------------------

.. automodule:: fe.stepactions.initializematerial
    :members: __doc__

``nodeforces`` - Concentrated node forces
-----------------------------------------

.. automodule:: fe.stepactions.nodeforces
    :members: __doc__

.. pprint:: fe.stepactions.nodeforces.documentation

``setfield`` - Set field to prescribed value
--------------------------------------------

.. automodule:: fe.stepactions.setfield
    :members: __doc__

.. pprint:: fe.stepactions.setfield.documentation

.. Solvers
.. *******

.. Constraints
.. ***********

.. Outputmanagers
.. **************

