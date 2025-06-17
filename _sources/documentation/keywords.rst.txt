.. list-table:: ``*advancedmaterial`` : definition of an advanced material
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - material description, multiline possible
    * - ``id``
      - ``string``
      - name of the property
    * - ``name``
      - ``string``
      - name of the property

.. list-table:: ``*analyticalfield`` : define an analytical field
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - definition
    * - ``name``
      - ``string``
      - name of analytical field
    * - ``type``
      - ``string``
      - type of analytical field (currently 'expression' only)

.. list-table:: ``*configureplots`` : customize the figures and axes
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - key=value pairs for configuration of figures and axes

.. list-table:: ``*constraint`` : define a constraint
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - definition of the constraint
    * - ``name``
      - ``string``
      - (optional) name of the constraint
    * - ``type``
      - ``string``
      - constraint type

.. list-table:: ``*element`` : definition of element(s)
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``numpy integer array``
      - Abaqus like element definiton lines
    * - ``elset``
      - ``string``
      - name
    * - ``provider``
      - ``string``
      - provider (library) for the element type. Default: Marmot
    * - ``type``
      - ``string``
      - assign one of the types definied in the elementlibrary

.. list-table:: ``*elset`` : definition of an element set
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - Abaqus like element set definiton lines
    * - ``elset``
      - ``string``
      - name
    * - ``generate``
      - ``string``
      - set True to generate from data line 1: start-element, end-element, step

.. list-table:: ``*exportplots`` : export your figures
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - key=value pairs for exporting of figures and axes

.. list-table:: ``*fieldoutput`` : define fieldoutput, which is used by outputmanagers
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - defintions lines for the output module
    * - ``jobname``
      - ``string``
      - (optional), name of job, standard=defaultJob

.. list-table:: ``*include`` : (optional) load extra .inp file (fragment), use relative path to current .inp
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``input``
      - ``string``
      - filename

.. list-table:: ``*job`` : definition of an analysis job
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``domain``
      - ``string``
      - define spatial domain: 1d, 2d, 3d
    * - ``name``
      - ``string``
      - (optional) name of job, standard = defaultJob
    * - ``starttime``
      - ``float``
      - (optional) start time of job, standard = 0.0

.. list-table:: ``*material`` : definition of a material
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``numpy float array``
      - material properties, multiline possible
    * - ``id``
      - ``string``
      - name of the property
    * - ``name``
      - ``string``
      - name of the property
    * - ``statevars``
      - ``integer``
      - (deprecated and ignored) number of statevars

.. list-table:: ``*modelgenerator`` : define a model generator, loaded from a module
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - key=value pairs
    * - ``generator``
      - ``string``
      - generator module
    * - ``name``
      - ``string``
      - (optional) name of the generator

.. list-table:: ``*node`` : definition of nodes
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``numpy float array``
      - Abaqus like node definiton lines: label, x, [y], [z]
    * - ``nset``
      - ``string``
      - name

.. list-table:: ``*nset`` : definition of an element set
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - Abaqus like node set definiton lines
    * - ``generate``
      - ``string``
      - set True to generate from data line 1: start-node, end-node, step
    * - ``nset``
      - ``string``
      - name

.. list-table:: ``*output`` : define an output module
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - defintions lines for the output module
    * - ``jobname``
      - ``string``
      - (optional), name of job, standard=defaultJob
    * - ``name``
      - ``string``
      - (optional), name of manager, standard=None
    * - ``type``
      - ``string``
      - output module

.. list-table:: ``*parameteridentification`` : identify material parameter for given x and y data
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - key=value pairs
    * - ``file``
      - ``string``
      - filename where identified parameters are written
    * - ``id``
      - ``string``
      - name of the property
    * - ``jobname``
      - ``string``
      - name of job
    * - ``plot``
      - ``string``
      - True|False plot result with final parameters
    * - ``xdata``
      - ``string``
      - filename where xData is given
    * - ``ydata``
      - ``string``
      - filename where yData is given

.. list-table:: ``*section`` : definition of an section
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - list of associated element sets
    * - ``material``
      - ``string``
      - associated id of defined material
    * - ``name``
      - ``string``
      - name
    * - ``thickness``
      - ``float``
      - associated element set
    * - ``type``
      - ``string``
      - type of the section

.. list-table:: ``*solver`` : definition of a solver
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``name``
      - ``string``
      - name of this solver
    * - ``solver``
      - ``string``
      - solver type
    * - ``data``
      - ``string``
      - (optional) define options which are passed to the respective solver instance.

.. list-table:: ``*step`` : definition of job steps
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``criticaliter``
      - ``integer``
      - maximum number of iterations to prevent from increasing the increment
    * - ``data``
      - ``string``
      - define step actions, which are handled by the corresponding stepaction modules
    * - ``jobname``
      - ``string``
      - (optional), name of job, standard=defaultJob
    * - ``maxinc``
      - ``float``
      - maximum size of increment
    * - ``maxiter``
      - ``integer``
      - maximum number of iterations
    * - ``maxnuminc``
      - ``integer``
      - maximum number of increments
    * - ``mininc``
      - ``float``
      - minimum size of increment
    * - ``steplength``
      - ``float``
      - time period of step

.. list-table:: ``*surface`` : definition of surface set
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - Abaqus like definition. Type 'element': elSet, faceID
    * - ``name``
      - ``string``
      - name
    * - ``type``
      - ``string``
      - type of surface (currently 'element' only)

.. list-table:: ``*updateconfiguration`` : update an configuration
    :width: 100%
    :widths: 25 25 40
    :header-rows: 1

    * - Option
      - Type
      - Description
    * - ``configuration``
      - ``string``
      -  name of the modified settings category
    * - ``data``
      - ``string``
      - key=value pairs

