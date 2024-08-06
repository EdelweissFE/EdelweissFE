Elements
========

Relevant module: ``edelweissfe.config.elementlibrary``

.. automodule:: edelweissfe.config.elementlibrary
   :members:

Provider ``marmot``
-------------------

Relevant module: ``edelweissfe.elements.marmotelement.element``

.. autoclass:: edelweissfe.elements.marmotelement.element.MarmotElementWrapper
   :members:

Provider ``marmotsingleqpelement``
----------------------------------

Relevant module: ``edelweissfe.elements.marmotsingleqpelement.element``

.. autoclass:: edelweissfe.elements.marmotsingleqpelement.element.MarmotMaterialWrappingElement
   :members:

.. literalinclude:: ../../../testfiles/QPMarmotMaterialHypoElastic/test.inp
    :language: edelweiss
    :caption: Example: ``testfiles/QPMarmotMaterialHypoElastic/test.inp``

Provider ``displacementelement``
--------------------------------

Relevant module: ``edelweissfe.elements.displacementelement.element``

.. autoclass:: edelweissfe.elements.displacementelement.element.DisplacementElement
   :members:

The elementType definition has to be of the following form:

**Elements**

One of the following element types needs to be included in the definition (``elementType[0:5]``).

- **Quad4**    - quadrilateral 2D element with 4 nodes.
- **Quad8**    - quadrilateral 2D element with 8 nodes.
- **Hexa8**    - hexahedron 3D element with 8 nodes.

**additional Parameters**

The following optional Parameters are also included in the element type definition.

- **R**     - reduced integration for element, in ``elementType[5]``, currently possible for Q8 and Q4.
- **E**     - extended integration for element, in ``elementType[5]``, currently possible for Q4.
- **PE**    - use plane strain for 2D elements, in ``elementType[6:8]`` or ``elementType[5:7]``.
- **PS**    - use plane stress for 2D elements, in ``elementType[6:8]`` or ``elementType[5:7]``.

.. note::
    If R or E is not given by the user, we assume regular integration.

    If PE or PS is not given by the user, we assume PE.

Example of the usage of a quadrilateral element with 4 nodes, regular integration and plane strain.

.. literalinclude:: ../../../examples/CantileverBeam/test4.inp
    :language: edelweiss
    :caption: Example: ``examples/CantileverBeam/test4.inp``

Example of the usage of a quadrilateral element with 8 nodes, regular integration and plane strain for the same example.

.. literalinclude:: ../../../examples/CantileverBeam/test8.inp
    :language: edelweiss
    :caption: Example: ``examples/CantileverBeam/test8.inp``

Example of the usage of a 3D hexahedron element with 8 nodes.

.. literalinclude:: ../../../examples/CantileverBeam/test3D.inp
    :language: edelweiss
    :caption: Example: ``examples/CantileverBeam/test3D.inp``

Implementing your own elements
------------------------------

Relevant module: ``edelweissfe.elements.base.baseelement``

.. automodule:: edelweissfe.elements.base.baseelement
   :members:
