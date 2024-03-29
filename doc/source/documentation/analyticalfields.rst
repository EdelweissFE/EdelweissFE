Analytical fields
=================

.. automodule:: fe.config.analyticalfields
    :members: __doc__

``scalarexpression`` - Field defined by expression
--------------------------------------------------

Relevant module ``fe.analyticalfields.scalarexpression``

.. automodule:: fe.analyticalfields.scalarexpression
    :members: __doc__

.. pprint:: fe.analyticalfields.scalarexpression.documentation
    :caption: Options

.. literalinclude:: ../../../testfiles/AnalyticalFieldsScalarExpression/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/AnalyticalFieldsScalarExpression/test.inp``

``randomscalar`` - A random field
---------------------------------

Relevant module ``fe.analyticalfields.randomscalar``

.. automodule:: fe.analyticalfields.randomscalar
    :members: __doc__

.. pprint:: fe.analyticalfields.randomscalar.documentation
    :caption: Options

.. literalinclude:: ../../../testfiles/AnalyticalFieldsRandomScalar/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/AnalyticalFieldsRandomScalar/test.inp``

Implementing your own fields
----------------------------

Subclass from the field base class in module ``fe.analyticalfields.base.analyticalfieldbase``

.. automodule:: fe.analyticalfields.base.analyticalfieldbase    
   :members: 
