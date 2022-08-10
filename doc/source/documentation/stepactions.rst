Step actions
============

.. automodule:: fe.config.stepactions
   :members: __doc__

``bodyforce`` - Bodyforce loads 
-------------------------------

Relevant module: ``fe.stepactions.bodyforce``

.. automodule:: fe.stepactions.bodyforce
   :members: __doc__

.. pprint:: fe.stepactions.bodyforce.documentation
   :caption: Options:

.. literalinclude:: ../../../testfiles/BodyForce3D/test.inp
   :language: console 
   :caption: Example: ``testfiles/BodyForce3D/test.inp``

``dirichlet`` - Standard birichlet BC
-------------------------------------

Relevant module: ``fe.stepactions.dirichlet``

.. automodule:: fe.stepactions.dirichlet
   :members: __doc__

.. pprint:: fe.stepactions.dirichlet.documentation
   :caption: Options:

``distributedload`` - surfaces loads
------------------------------------

Relevant module: ``fe.stepactions.distributedload``

.. automodule:: fe.stepactions.distributedload
    :members: __doc__

.. pprint:: fe.stepactions.distributedload.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/DLoad/test.inp
   :language: console 
   :caption: Example: ``testfiles/DLoad/test.inp``

``geostatic`` - Geostatic stress definition
-------------------------------------------

Relevant module: ``fe.stepactions.geostatic``

.. automodule:: fe.stepactions.geostatic
    :members: __doc__

.. pprint:: fe.stepactions.geostatic.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/GeoStatic/test.inp
   :language: console 
   :caption: Example: ``testfiles/GeoStatic/test.inp``

``indirectcontrol`` - Indirect displacement control
---------------------------------------------------

Relevant module: ``fe.stepactions.indirectcontrol``

.. automodule:: fe.stepactions.indirectcontrol
    :members: __doc__

.. pprint:: fe.stepactions.indirectcontrol.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/IndirectDisplacementControl//test.inp
   :language: console
   :caption: Example: ``testfiles/IndirectDisplacementControl/test.inp``

``indirectcontractioncontrol`` - Indirect displacement -- contraction ring control
----------------------------------------------------------------------------------

Relevant module: ``fe.stepactions.indirectcontrol``

.. automodule:: fe.stepactions.indirectcontractioncontrol
    :members: __doc__

.. pprint:: fe.stepactions.indirectcontractioncontrol.documentation
   :caption: Options:

``initalizemateiral`` - Initialize materials
--------------------------------------------

Relevant module: ``fe.stepactions.initializematerial``

.. automodule:: fe.stepactions.initializematerial
    :members: __doc__

``nodeforces`` - Concentrated node forces
-----------------------------------------

Relevant module: ``fe.stepactions.nodeforces``

.. automodule:: fe.stepactions.nodeforces
    :members: __doc__

.. pprint:: fe.stepactions.nodeforces.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/NodeForces/test.inp
   :language: console 
   :caption: Example: ``testfiles/NodeForces/test.inp``

``setfield`` - Set field to prescribed value
--------------------------------------------

Relevant module: ``fe.stepactions.setfield``

.. automodule:: fe.stepactions.setfield
    :members: __doc__

.. pprint:: fe.stepactions.setfield.documentation
   :caption: Options:


``modelupdate`` - Update the model
----------------------------------

Relevant module: ``fe.stepactions.modelupdate``

.. automodule:: fe.stepactions.modelupdate
    :members: __doc__

.. pprint:: fe.stepactions.modelupdate.documentation
   :caption: Options:

Implementing your own step actions
----------------------------------

Subclass from the step action base class in module ``fe.stepactions.base.stepactionbase``

.. automodule:: fe.stepactions.base.stepactionbase    
   :members: 
