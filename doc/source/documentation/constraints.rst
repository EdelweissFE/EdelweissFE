Constraints
===========

.. automodule:: fe.config.constraints
    :members: __doc__

``equalvaluelagrangian`` - Constrain nodal values to equal values 
-----------------------------------------------------------------

Module ``fe.constraints.equalvaluelagrangian``

.. automodule:: fe.constraints.equalvaluelagrangian
    :members: __doc__

.. pprint:: fe.constraints.equalvaluelagrangian.documentation
    :caption: Options:

.. literalinclude:: ../../../testfiles/EqualValueLagrangianConstraint/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/EqualValueLagrangianConstraint/test.inp``


``equalvaluepenalty`` - Constrain nodal values to equal values 
--------------------------------------------------------------

Module ``fe.constraints.equalvaluepenalty``

.. automodule:: fe.constraints.equalvaluepenalty
    :members: __doc__

.. pprint:: fe.constraints.equalvaluepenalty.documentation
    :caption: Options:

.. literalinclude:: ../../../testfiles/EqualValuePenaltyConstraint/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/EqualValuePenaltyConstraint/test.inp``


``linearizedrigidbody`` - Linearized rigid body constraints in 2D
------------------------------------------------------------------

Module ``fe.constraints.linearizedrigidbody``

.. automodule:: fe.constraints.linearizedrigidbody
    :members: __doc__

.. pprint:: fe.constraints.linearizedrigidbody.documentation
    :caption: Options:

.. literalinclude:: ../../../testfiles/LinearizedRigidBodyConstraint/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/LinearizedRigidBodyConstraint/test.inp``


``rigidbody`` - Geometrically exact rigid body constraints in 3D
---------------------------------------------------------------------------

Module ``fe.constraints.rigidbody``

.. automodule:: fe.constraints.rigidbody
    :members: __doc__

.. pprint:: fe.constraints.rigidbody.documentation
    :caption: Options

.. literalinclude:: ../../../testfiles/RigidBodyConstraintLargeDeformations3D/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/RigidBodyConstraintLargeDeformations3D/test.inp``

``penaltyindirectcontrol`` - Penalty based indirect control
-----------------------------------------------------------

Module ``fe.constraints.penaltyindirectcontrol``

.. automodule:: fe.constraints.penaltyindirectcontrol
    :members: __doc__

.. pprint:: fe.constraints.penaltyindirectcontrol.documentation
    :caption: Options

.. literalinclude:: ../../../testfiles/PenaltyBasedIndirectControl/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/PenaltyBasedIndirectControl/test.inp``


Implementing your own constraints
---------------------------------

Subclass from the constraint base class in module ``fe.constraints.base.constraintbase``

.. automodule:: fe.constraints.base.constraintbase
    :members: 
