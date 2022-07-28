``*analyticalfield`` : define an analytical field
 
.. list-table:: Options
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
 
``*configureplots`` : customize the figures and axes
 
.. list-table:: Options
    :widths: 25 25 40
    :header-rows: 1
 
    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - key=value pairs for configuration of figures and axes
 
``*constraint`` : define a constraint
 
.. list-table:: Options
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
 
``*element`` : definition of element(s)
 
.. list-table:: Options
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
 
``*elset`` : definition of an element set
 
.. list-table:: Options
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
 
``*exportplots`` : export your figures
 
.. list-table:: Options
    :widths: 25 25 40
    :header-rows: 1
 
    * - Option
      - Type
      - Description
    * - ``data``
      - ``string``
      - key=value pairs for exporting of figures and axes
 
``*fieldoutput`` : define fieldoutput, which is used by outputmanagers
 
.. list-table:: Options
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
 
``*include`` : (optional) load extra .inp file (fragment), use relative path to current .inp
 
.. list-table:: Options
    :widths: 25 25 40
    :header-rows: 1
 
    * - Option
      - Type
      - Description
    * - ``input``
      - ``string``
      - filename
 
``*job`` : definition of an analysis job
 
.. list-table:: Options
    :widths: 25 25 40
    :header-rows: 1
 
    * - Option
      - Type
      - Description
    * - ``domain``
      - ``string``
      - define spatial domain: 1d, 2d, 3d
    * - ``linsolver``
      - ``string``
      - (optional) define linear solver, standard = superlu
    * - ``name``
      - ``string``
      - (optional) name of job, standard = defaultJob
    * - ``solver``
      - ``string``
      - (optional) define solver, standard = NIST
    * - ``starttime``
      - ``float``
      - (optional) start time of job, standard = 0.0
 
``*material`` : definition of a material
 
.. list-table:: Options
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
 
``*modelgenerator`` : define a model generator, loaded from a module
 
.. list-table:: Options
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
 
``*node`` : definition of nodes
 
.. list-table:: Options
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
 
``*nset`` : definition of an element set
 
.. list-table:: Options
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
 
``*output`` : define an output module
 
.. list-table:: Options
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
 
``*parameteridentification`` : identify material parameter for given x and y data
 
.. list-table:: Options
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
 
``*section`` : definition of an section
 
.. list-table:: Options
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
 
``*step`` : definition of job steps
 
.. list-table:: Options
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
 
``*surface`` : definition of surface set
 
.. list-table:: Options
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
 
``*updateconfiguration`` : update an configuration
 
.. list-table:: Options
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
 
