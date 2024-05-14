Materials
=========

Materials are defined in EdelweissFE using the ``*material`` keyword.
Mandatory arguments are

.. list-table:: ``*material`` : definition of a material
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``id``
      - ``string``
      - A unique ID, which is used for referencing the material in EdelweissFE.
    * - ``name``
      - ``string``
      - The name of of the material.
    * - ``datalines``
      - ``numpy float array``
      - The material properties as a float vector, multiline possible.


Materials are assigned to elements by means of :ref:`sections`.
