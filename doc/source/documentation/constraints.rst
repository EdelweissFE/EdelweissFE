Constraints
===========

.. automodule:: edelweissfe.config.constraints
    :members: __doc__

``equalvaluelagrangian`` - Constrain nodal values to equal values
-----------------------------------------------------------------

Module ``edelweissfe.constraints.equalvaluelagrangian``

.. automodule:: edelweissfe.constraints.equalvaluelagrangian
    :members: __doc__

.. pprint:: edelweissfe.constraints.equalvaluelagrangian.documentation
    :caption: Options:

.. literalinclude:: ../../../testfiles/EqualValueLagrangianConstraint/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/EqualValueLagrangianConstraint/test.inp``


``equalvaluepenalty`` - Constrain nodal values to equal values
--------------------------------------------------------------

Module ``edelweissfe.constraints.equalvaluepenalty``

.. automodule:: edelweissfe.constraints.equalvaluepenalty
    :members: __doc__

.. pprint:: edelweissfe.constraints.equalvaluepenalty.documentation
    :caption: Options:

.. literalinclude:: ../../../testfiles/EqualValuePenaltyConstraint/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/EqualValuePenaltyConstraint/test.inp``


``linearizedrigidbody`` - Linearized rigid body constraints
-----------------------------------------------------------

Module ``edelweissfe.constraints.linearizedrigidbody``

.. automodule:: edelweissfe.constraints.linearizedrigidbody
    :members: __doc__

.. pprint:: edelweissfe.constraints.linearizedrigidbody.documentation
    :caption: Options:

.. literalinclude:: ../../../testfiles/LinearizedRigidBodyConstraint/test.inp
    :language: edelweiss
    :caption: Example 2D: ``testfiles/LinearizedRigidBodyConstraint/test.inp``

.. literalinclude:: ../../../testfiles/LinearizedRigidBodyConstraint2D/test.inp
    :language: edelweiss
    :caption: Example 2D: ``testfiles/LinearizedRigidBodyConstraint2D/test.inp``

.. literalinclude:: ../../../testfiles/LinearizedRigidBodyConstraint3D/test.inp
    :language: edelweiss
    :caption: Example 3D: ``testfiles/LinearizedRigidBodyConstraint3D/test.inp``


``rigidbody`` - Geometrically exact rigid body constraints in 3D
---------------------------------------------------------------------------

Module ``edelweissfe.constraints.rigidbody``

.. automodule:: edelweissfe.constraints.rigidbody
    :members: __doc__

.. pprint:: edelweissfe.constraints.rigidbody.documentation
    :caption: Options

.. literalinclude:: ../../../testfiles/RigidBodyConstraintLargeDeformations3D/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/RigidBodyConstraintLargeDeformations3D/test.inp``

``penaltyindirectcontrol`` - Penalty based indirect control
-----------------------------------------------------------

Module ``edelweissfe.constraints.penaltyindirectcontrol``

.. automodule:: edelweissfe.constraints.penaltyindirectcontrol
    :members: __doc__

.. pprint:: edelweissfe.constraints.penaltyindirectcontrol.documentation
    :caption: Options

.. literalinclude:: ../../../testfiles/PenaltyBasedIndirectControl/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/PenaltyBasedIndirectControl/test.inp``


Implementing your own constraints
---------------------------------

Subclass from the constraint base class in module ``edelweissfe.constraints.base.constraintbase``

.. automodule:: edelweissfe.constraints.base.constraintbase
    :members:
