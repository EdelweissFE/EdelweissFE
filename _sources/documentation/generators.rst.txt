Meshgenerators
==============

Generate meshes on the fly. Usage:

.. code-block:: console

    *job, name=job, domain=3d, solver=NISTParallel
    
    *modelGenerator, generator=boxGen, name=gen
    nX      =4
    nY      =8
    nZ      =2
    lX      =20
    lY      =40
    lZ      =1
    elType  =C3D20R


``boxgen`` - A 3D box mesh generator
------------------------------------

.. automodule:: fe.generators.boxgen
   :members: __doc__

.. pprint:: fe.generators.boxgen.documentation
   
Example: ``testfiles/BoxGen/test.inp``

.. literalinclude:: ../../../testfiles/BoxGen/test.inp
   :language: console

``planerectquad`` - A 2D plane rectangular mesh generator
---------------------------------------------------------

.. automodule:: fe.generators.planerectquad
   :members: __doc__

.. pprint:: fe.generators.planerectquad.documentation

Example: ``testfiles/NodeForces/test.inp``

.. literalinclude:: ../../../testfiles/NodeForces/test.inp
   :language: console

``cubit`` - A cubit mesh generator
---------------------------------------------------------

.. automodule:: fe.generators.cubit
   :members: __doc__

.. pprint:: fe.generators.cubit.documentation

Example: ``testfiles/CubitGen/test.inp``

.. literalinclude:: ../../../testfiles/CubitGen/test.inp
   :language: console
