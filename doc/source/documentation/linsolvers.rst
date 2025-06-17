Linear solvers
==============

Linear solvers are defined in EdelweissFE after the ``*solver`` keyword using ``linsolver`` and an optional configuration file ``linsolverConfigFile`` as a data line.
The ``linsolverConfigFile`` needs to be in ``.json`` format.

Choose a linsolver after the ``*solver`` keyword:

.. code-block:: edelweiss

    *solver, solver=NIST, name=theSolver
    linsolver=gmres
    linsolverConfigFile=opt.json

.. list-table:: Currently available linear solvers
    :width: 100%
    :widths: 15 1 25
    :header-rows: 1

    * - Name
      - Direct solver
      - Relevant module
    * - ``superlu``
      - ✓
      - ``scipy.sparse.linalg.spsolve``
    * - ``umfpack``
      - ✓
      - ``scipy.sparse.linalg.spsolve``
    * - ``pardiso``
      - ✓
      - ``edelweissfe.linsolve.pardiso.pardiso``
    * - ``panuapardiso``
      - ✓
      - ``edelweissfe.linsolve.panuapardiso.panuapardiso``
    * - ``klu``
      - ✓
      - ``edelweissfe.linsolve.klu.klu``
    * - ``petsclu``
      - ✓
      - ``edelweissfe.linsolve.petsclu.petsclu``
    * - ``mumps``
      - ✓
      - ``edelweissfe.linsolve.mumps.mumps``
    * - ``gmres``
      - ✗
      - ``edelweissfe.linsolve.gmres.gmres``

Currently, only the linsolver ``gmres`` allows the use of an optional configuration file ``linsolverConfigFile``.

Choose the options for the linsolver (in this case ``gmres``) in an extra file:

.. code-block:: json

    	{
	"precondopts":
	{
	"presmoother": ["block_gauss_seidel", {"iterations": 15}],
	"postsmoother": ["block_gauss_seidel", {"iterations": 15}],
	},
	"linsolveopts": {"maxiter": 1, "restart": 1500}
	}
