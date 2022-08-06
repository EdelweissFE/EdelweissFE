Constraints
===========

.. automodule:: fe.config.constraints
    :members: __doc__

``linearizedrigidbody`` - Linearized rigid body constraints in 2D
------------------------------------------------------------------

Module ``fe.constraints.linearizedrigidbody``

.. automodule:: fe.constraints.linearizedrigidbody
    :members: __doc__

.. pprint:: fe.constraints.linearizedrigidbody.documentation
    :caption: Options:

.. literalinclude:: ../../../testfiles/LinearizedRigidBodyConstraint/test.inp
    :language: console
    :caption: Example: ``testfiles/LinearizedRigidBodyConstraint/test.inp``


``rigidbody`` - Geometrically exact rigid body constraints in 3D
---------------------------------------------------------------------------

Module ``fe.constraints.rigidbody``

.. automodule:: fe.constraints.rigidbody
    :members: __doc__

.. pprint:: fe.constraints.rigidbody.documentation
    :caption: Options

.. literalinclude:: ../../../testfiles/RigidBodyConstraintLargeDeformations3D/test.inp
    :language: console
    :caption: Example: ``testfiles/RigidBodyConstraintLargeDeformations3D/test.inp``


Implementing your own constraints
---------------------------------

Subclass from the constraint base class in module ``fe.constraints.base.constraintbase``

.. automodule:: fe.constraints.base.constraintbase
    :members: 
