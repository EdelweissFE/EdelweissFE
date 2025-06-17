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

Provider ``edelweiss``
----------------------

The element type definition has to be of the following form:

**Elements**

One of the following element types needs to be included in the definition (``type`` in .inp file).

- **CPE4**    - quadrilateral 2D element with 4 nodes and plane strain.
- **CPE8**    - quadrilateral 2D element with 8 nodes and plane strain.
- **CPS4**    - quadrilateral 2D element with 4 nodes and plane stress.
- **CPS8**    - quadrilateral 2D element with 8 nodes and plane stress.
- **C3D8**    - hexahedron 3D element with 8 nodes.
- **C3D20**    - hexahedron 3D element with 20 nodes.

**additional Parameters**

The following optional Parameters are also included in the element type definition:

- **R**     - reduced integration for element, at the end of elementType.
- **E**     - extended integration for element, at the end of elementType.
- **N**     - (optional) regular integration, at the end of elementType.
- **TL**    - use the total Lagrangian element, before integration type.

Geometrically linear element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This element uses a geometrically linear approach for the computation.

Relevant module: ``edelweissfe.elements.displacementelement.element``

.. autoclass:: edelweissfe.elements.displacementelement.element.DisplacementElement
   :members:

.. note::
    If R or E is not given by the user, we assume regular integration.

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

.. note::
    This element is only compatible with the material provider ``edelweiss``!

.. note::
    This element works with the Cauchy stress :math:`\sigma` and the tangent :math:`\frac{d\sigma}{d\varepsilon}`, where
    :math:`\varepsilon` stands for the linearized strain tensor!

Total Lagrange element
~~~~~~~~~~~~~~~~~~~~~~

This element is geometrically nonlinear and uses the total Lagrange algorithm for finite strain.

Relevant module: ``edelweissfe.elements.displacementtlelement.element``

.. autoclass:: edelweissfe.elements.displacementtlelement.element.DisplacementTLElement
    :members:

.. note::
    This element is only compatible with the material provider ``edelweiss``!

.. note::
    For hyperelastic materials, this element works with the second Piola-Kirchhoff stress :math:`\mathbf{S}` and the tangent
    :math:`\frac{d\mathbf{S}}{d\mathbf{E}}`, where :math:`\mathbf{E}` stands for the Green-Lagrange strain tensor!

    For other materials, this element works with the Kirchhoff stress :math:`\tau` and the tangent
    :math:`\frac{d\tau}{d\mathbf{F}}`, where :math:`\mathbf{F}` stands for the deformation gradient!

Implementing your own elements
------------------------------

Relevant module: ``edelweissfe.elements.base.baseelement``

.. automodule:: edelweissfe.elements.base.baseelement
   :members: