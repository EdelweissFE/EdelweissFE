Parallelization
===============

EdelweissFE makes use of OpenMP to parallelize the computation of finite elements and for certain solvers, such as SuperLU, UMFPACK, or PARDISO.

If a parallel solver (e.g, NISTParallel, NISTPArcLength) is selected in the .inp file, EdelweissFE  automatically determines the maximum number of threads,
depending on the host architecture.
However, it is RECOMMENDED to enforce a fixed number of threads by running

.. code-block:: console

    OMP_NUM_THREADS=XX python edelweiss.py INPUT.inp

This ensures that the same number of threads ``XX`` is employed both in EdelweissFE as well as in the underlying Intel MKL (e.g., if the PARDISO linear solver is used).
