Constraints
===========

Constraints are defined globally,
and may be used for introducing additional scalar
equations into the global equation system.

.. code-block:: console

    *constraint, type=rigidbody, name=rb1
    nSet=right
    referencePoint=rBottom

``linearlizedrigidbody`` - Linearized rigid body constraints in 2D
------------------------------------------------------------------

.. automodule:: fe.constraints.linearizedrigidbody
   :members: __doc__

.. pprint:: fe.constraints.linearizedrigidbody.documentation

``rigidbody`` - Geometrically exact rigid body constraints in 3D
---------------------------------------------------------------------------

.. automodule:: fe.constraints.rigidbody
   :members: __doc__

.. pprint:: fe.constraints.rigidbody.documentation
