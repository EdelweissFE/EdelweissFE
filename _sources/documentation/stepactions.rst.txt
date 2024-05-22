Step actions
============

.. automodule:: edelweissfe.config.stepactions
   :members: __doc__

``bodyforce`` - Bodyforce loads
-------------------------------

Relevant module: ``edelweissfe.stepactions.bodyforce``

.. automodule:: edelweissfe.stepactions.bodyforce
   :members: __doc__

.. pprint:: edelweissfe.stepactions.bodyforce.documentation
   :caption: Options:

.. literalinclude:: ../../../testfiles/BodyForce3D/test.inp
   :language: edelweiss
   :caption: Example: ``testfiles/BodyForce3D/test.inp``

``dirichlet`` - Standard dirichlet BC
-------------------------------------

Relevant module: ``edelweissfe.stepactions.dirichlet``

.. automodule:: edelweissfe.stepactions.dirichlet
   :members: __doc__

.. pprint:: edelweissfe.stepactions.dirichlet.documentation
   :caption: Options:

``distributedload`` - Distributed surface loads
-----------------------------------------------

Relevant module: ``edelweissfe.stepactions.distributedload``

.. automodule:: edelweissfe.stepactions.distributedload
    :members: __doc__

.. pprint:: edelweissfe.stepactions.distributedload.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/DLoad/test.inp
   :language: edelweiss
   :caption: Example: ``testfiles/DLoad/test.inp``

``geostatic`` - Geostatic stress states
---------------------------------------

Relevant module: ``edelweissfe.stepactions.geostatic``

.. automodule:: edelweissfe.stepactions.geostatic
    :members: __doc__

.. pprint:: edelweissfe.stepactions.geostatic.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/GeoStatic/test.inp
   :language: edelweiss
   :caption: Example: ``testfiles/GeoStatic/test.inp``

``indirectcontrol`` - Indirect displacement control
---------------------------------------------------

Relevant module: ``edelweissfe.stepactions.indirectcontrol``

.. automodule:: edelweissfe.stepactions.indirectcontrol
    :members: __doc__

.. pprint:: edelweissfe.stepactions.indirectcontrol.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/IndirectDisplacementControl//test.inp
   :language: edelweiss
   :caption: Example: ``testfiles/IndirectDisplacementControl/test.inp``

``indirectcontractioncontrol`` - Indirect displacement -- contraction ring control
----------------------------------------------------------------------------------

Relevant module: ``edelweissfe.stepactions.indirectcontrol``

.. automodule:: edelweissfe.stepactions.indirectcontractioncontrol
    :members: __doc__

.. pprint:: edelweissfe.stepactions.indirectcontractioncontrol.documentation
   :caption: Options:

``initializematerial`` - Initialize materials
---------------------------------------------

Relevant module: ``edelweissfe.stepactions.initializematerial``

.. automodule:: edelweissfe.stepactions.initializematerial
    :members: __doc__

``modelupdate`` - Update the model
----------------------------------

Relevant module: ``edelweissfe.stepactions.modelupdate``

.. automodule:: edelweissfe.stepactions.modelupdate
    :members: __doc__

.. pprint:: edelweissfe.stepactions.modelupdate.documentation
   :caption: Options:

``nodeforces`` - Concentrated node forces
-----------------------------------------

Relevant module: ``edelweissfe.stepactions.nodeforces``

.. automodule:: edelweissfe.stepactions.nodeforces
    :members: __doc__

.. pprint:: edelweissfe.stepactions.nodeforces.documentation
   :caption: Options:


.. literalinclude:: ../../../testfiles/NodeForces/test.inp
   :language: edelweiss
   :caption: Example: ``testfiles/NodeForces/test.inp``

``setfield`` - Set a field to a prescribed value
------------------------------------------------

Relevant module: ``edelweissfe.stepactions.setfield``

.. automodule:: edelweissfe.stepactions.setfield
    :members: __doc__

.. pprint:: edelweissfe.stepactions.setfield.documentation
   :caption: Options:

``setinitialconditions`` - Set initial conditions to elements
-------------------------------------------------------------

Relevant module: ``edelweissfe.stepactions.setinitialconditions``

.. automodule:: edelweissfe.stepactions.setinitialconditions
    :members: __doc__

.. pprint:: edelweissfe.stepactions.setinitialconditions.documentation
   :caption: Options:

Implementing your own step actions
----------------------------------

Subclass from the step action base class in module ``edelweissfe.stepactions.base.stepactionbase``

.. automodule:: edelweissfe.stepactions.base.stepactionbase
   :members:
