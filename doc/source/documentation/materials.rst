Materials
=========

Relevant module: ``edelweissfe.config.materiallibrary``

Materials are defined in EdelweissFE using the ``*material`` or ``*advancedmaterial`` keyword.
Arguments for the ``*material`` keyword are

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
    * - ``provider``
      - ``string``
      - (optional) The material provider.

The material properties for materials using the ``*material`` keyword are assigned as seen here:

.. code-block:: edelweiss
    :caption: Use different non-advanced materials. Example:

    *material, name=neohookewaplastic, id=Mat1, provider=edelweissmaterial
    91304.34783, 100000., 260, 70, 320, 9

Arguments for the ``*advancedmaterial`` keyword are

.. list-table:: ``*advancedmaterial`` : definition of a material
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
      - ``string``
      - The material properties as a vector of strings, multiline possible.
    * - ``provider``
      - ``string``
      - (optional) The material provider.

The material properties for materials using the ``*advancedmaterial`` keyword are assigned as seen here:

.. code-block:: edelweiss
    :caption: Use different advanced materials. Example:

    *advancedmaterial, name=hyperplasticadvanced, id=hyperplasticadvanced, provider=edelweissmaterial
    psi_e='mu/2 * (I1/J**(2/3) - 3) + K/8 * (J**2 + 1/J**2 - 2)'
    mu=91304.34783, K=100000., fy0=260, HLin=70, dfy=320, delta=9
    a='2,3,4'

If the material provider is not given, ``marmotmaterial`` is assumed. Materials are assigned to elements by means of :ref:`sections`.

Provider ``marmotmaterial``
---------------------------

Relevant module: `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_.

Provider ``edelweiss``
----------------------

Relevant module: ``edelweissfe.materials``

.. list-table:: Currently available materials with this provider
    :width: 100%
    :widths: 15 25 15
    :header-rows: 1

    * - Name
      - Description
      - Keyword
    * - ``linearelastic``
      - Linear elastic material.
      - ``*material``
    * - ``vonmises``
      - Von Mises material.
      - ``*material``
    * - ``neohookewa``
      - Neo-Hookean Pence-Gou formulation 'a' material.
      - ``*material``
    * - ``neohookewb``
      - Neo-Hookean Pence-Gou formulation 'b' material.
      - ``*material``
    * - ``neohookewc``
      - Neo-Hookean Pence-Gou formulation 'c' material.
      - ``*material``
    * - ``hyperelasticadvanced``
      - Hyperelastic material with advanced defined energy density function using I1 and J.
      - ``*advancedmaterial``
    * - ``hyperelasticadvancedi2extended``
      - Hyperelastic material with advanced defined energy density function using I1, I2, J and C itself.
      - ``*advancedmaterial``
    * - ``neohookewaplastic``
      - Neo-Hookean Pence-Gou formulation 'a' material with J2 plasticity.
      - ``*material``
    * - ``neohookewbplastic``
      - Neo-Hookean Pence-Gou formulation 'b' material with J2 plasticity.
      - ``*material``
    * - ``neohookewcplastic``
      - Neo-Hookean Pence-Gou formulation 'c' material with J2 plasticity.
      - ``*material``
    * - ``hyperplasticadvanced``
      - Hyperelastic-plastic material with advanced defined energy density function using I1 and J.
      - ``*advancedmaterial``

Linear elastic material
~~~~~~~~~~~~~~~~~~~~~~~

This material uses a linear elastic relation between stress and strain.
The linear elastic material can be used as 2D plane stress and plane strain material as well as 3D material.
This material needs the following material properties in the correct order

#. **E**       - Elasticity module (Young's modulus).
#. :math:`\mathbf{\nu}`       - Poisson's ratio.

For the 2D plane strain and 3D material the following law in voigt notation is used

.. math::
   \begin{bmatrix}\sigma_{11}\\ \sigma_{22} \\ \sigma_{33}\\ \sigma_{12}\\ \sigma_{23}\\ \sigma_{13}\end{bmatrix} =
   \frac{E}{(1+\nu)(1-2\nu)}
   \begin{bmatrix}(1-\nu) & \nu & \nu & 0 & 0 & 0 \\
                  & (1-\nu) & \nu & 0 & 0 & 0\\
                  & & (1-\nu) & 0 & 0 & 0\\
                  & & & \frac{1-2\nu}{2} & 0 & 0\\
                  & & & & \frac{1-2\nu}{2} & 0\\
                  \text{symm.}& & & & & \frac{1-2\nu}{2}\end{bmatrix}
   \begin{bmatrix}\varepsilon_{11}\\ \varepsilon_{22} \\ \varepsilon_{33}\\ \gamma_{12}\\ \gamma_{23}\\ \gamma_{13}\end{bmatrix}

and for the 2D plane stress material the following law is used

.. math::
   \begin{bmatrix}\sigma_{11}\\ \sigma_{22} \\ \sigma_{12}\end{bmatrix} =
   \frac{E}{1-\nu^2}
   \begin{bmatrix}1 & \nu & 0\\
                  & 1 & 0\\
                  \text{symm.}& & \frac{1-\nu}{2}\end{bmatrix}
   \begin{bmatrix}\varepsilon_{11}\\ \varepsilon_{22} \\ \gamma_{12}\end{bmatrix}.

For the second case the third strain component gets calculated by using

.. math::
   \varepsilon_{33} = -\frac{\nu}{1-\nu} (\varepsilon_{11}+\varepsilon_{22}).

.. autoclass:: edelweissfe.materials.linearelastic.linearelastic.LinearElasticMaterial
   :members:

Von Mises material
~~~~~~~~~~~~~~~~~~

This material uses the same linear elastic law as the linear elastic material until plasticity is reached.
The von Mises material can be used as 2D plane strain and 3D material.
This material needs the following material properties in the correct order

#. **E**       - Elasticity module (Young's modulus).
#. :math:`\mathbf{\nu}`          - Poisson's ratio.
#. :math:`\mathbf{f_{y0}}`       - Yield stress.
#. :math:`\mathbf{H_{lin}}`      - Linear plastic hardening parameter.
#. :math:`\mathbf{\Delta f_y}`   - Multiplicator for nonlinear isotropic hardening.
#. :math:`\mathbf{\delta}`       - Exponent for nonlinear isotropic hardening.

The law used for nonlinear isotropic hardening in this material is :math:`f_y(\kappa)=f_{y0}+H_{lin}\kappa+\Delta f_{y} e^{-\delta\kappa}`
with :math:`\kappa` as the hardening parameter.
Plasticity is reached once the yield function f is :math:`f(\kappa)=||\mathbf{s}||-\sqrt{\frac{2}{3}}f_y(\kappa) > 0` with :math:`\mathbf{s}`
as the deviatoric stress.
With plasticity reached the material calculates a new :math:`\Delta\kappa` using Newton's method which is afterwards added to
the hardening parameter :math:`\kappa`. In the end the full material tangent and the back projected stress get calculated. For the back
projected stress

.. math::
   \mathbf{\sigma} = \mathbf{\sigma}^{trial} - 2G\sqrt{\frac{3}{2}}\Delta\kappa\frac{\mathbf{s}}{||\mathbf{s}||}

is used with

.. math::
   G = \frac{E}{2(1 + \nu)}

and :math:`\mathbf{s}` as the deviatoric stress.

.. autoclass:: edelweissfe.materials.vonmises.vonmises.VonMisesMaterial
   :members:

Elastic Neo-Hookean W(I1, J) Pence-Gou [a] materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These materials need the following parameters as input data in the correct order:

#. :math:`\mu`       - Shear modulus.
#. :math:`\kappa`    - Bulk modulus.

These materials base on a energy density function :math:`W(I_1, J)`, where :math:`I_1=\text{trace}(\mathbf{b})` is the first
and :math:`J=\det(\mathbf{F})` is the second invariant of the left Cauchy-Green tensor :math:`\mathbf{b}=\mathbf{FF}^\text{T}`
with the deformation gradient :math:`\mathbf{F}`, the Kirchhoff stress is then:

.. math::
   \tau = 2 \frac{\partial W}{\partial I_1} \mathbf{b} + J\frac{\partial W}{\partial J} \mathbf{I}.

The tangent modulus is then given by

.. math::
   C^{\tau \text{F}}_{ijkL} = 2b_{ij}\frac{\partial^2 W}{\partial I_1\partial F_{kL}} + 2\frac{\partial W}{\partial I_1}
   \frac{\partial b_{ij}}{\partial F_{kL}} +\frac{\partial}{\partial F_{kL}}\left[J\frac{\partial W}{\partial J}\right] \delta_{ij}.

The :math:`\mathbf{W_a}` **material** is defined by

.. math::
   W_a (I_1, J) = \frac{\mu}{2} (I_1 - 3) + \left (\frac{\kappa}{2} - \frac{\mu}{3} \right) (J - 1)^2 - \mu \ln (J).

The stress can then be calculated with

.. math::
   \tau = \mu \mathbf{b} + \left[\left(\kappa - \frac{2\mu}{3}\right) (J^2-J) - \mu\right] \mathbf{I}

and the tangent modulus is

.. math::
   C^{\tau \text{F}}_{ijkL} = \mu\frac{\partial b_{ij}}{\partial F_{kL}} + \left(\kappa - \frac{2\mu}{3}\right) (2J^2-J) \delta_{ij} (F^{-1})_{Lk}.

.. autoclass:: edelweissfe.materials.neohooke.neohookepencegouformulationa.NeoHookeanWaMaterial
   :members:

The :math:`\mathbf{W_b}` **material** is defined by

.. math::
   W_b (I_1, J) = \frac{\mu}{2} \left(\frac{I_1}{J^{2/3}} - 3\right) + \frac{\kappa}{8} \left (J^2 + \frac{1}{J^2} - 2 \right).

The stress can then be calculated with

.. math::
   \tau = \frac{\mu}{J^{2/3}} \mathbf{b} + \left[ \frac{\kappa}{4} \left( J^2 - \frac{1}{J^2} \right) -\frac{\mu I_1}{3 J^{2/3}}\right] \mathbf{I}

and the tangent modulus is

.. math::
   C^{\tau \text{F}}_{ijkL} = \frac{\mu}{J^{2/3}}\frac{\partial b_{ij}}{\partial F_{kL}} - \frac{2\mu}{3J^{2/3}} b_{ij} (F^{-1})_{Lk} -
   \frac{2\mu}{3J^{2/3}}\delta_{ij} F_{kL} + \left[ \frac{2\mu I_1}{9 J^{2/3}} + \frac{\kappa}{2} \left( J^2 + \frac{1}{J^2} \right) \right] \delta_{ij}\left(F^{-1}\right)_{Lk}.

.. autoclass:: edelweissfe.materials.neohooke.neohookepencegouformulationb.NeoHookeanWbMaterial
   :members:

The :math:`\mathbf{W_c}` **material** is defined by

.. math::
   W_c (I_1, J) = \frac{\mu}{2} \left(I_1 - 3\right) + \frac{3\mu^2}{3\kappa - 2\mu} \left(J^{\frac{2}{3} - \frac{\kappa}{\mu}} - 1\right).

The stress can then be calculated with

.. math::
   \tau = \mu \left(\mathbf{b}-J^{\frac{2}{3} - \frac{\kappa}{\mu}}\mathbf{I}\right)

and the tangent modulus is

.. math::
   C^{\tau \text{F}}_{ijkL} = \mu\frac{\partial b_{ij}}{\partial F_{kL}} + J^{\frac{2}{3}-\frac{\kappa}{\mu}}\left(\kappa - \frac{2\mu}{3}\right) \delta_{ij} (F^{-1})_{Lk}.

.. autoclass:: edelweissfe.materials.neohooke.neohookepencegouformulationc.NeoHookeanWcMaterial
   :members:

Advanced hyperelastic materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These materials need the following parameters as input data using keywords:

#. ``mu``       - Shear modulus.
#. ``K``        - Bulk modulus.
#. ``psi_e``    - Energy density function.
#. ``a``        - (Optional) your own parameters in ``psi_e`` (assignment with commas and quotation mark: ``a = '1.,2.,3'``).

The right Cauchy-Green tensor is defined as :math:`\mathbf{C}=\mathbf{F}^\text{T}\mathbf{F}`. These materials use the Kirchhoff stress and its corresponding tangent modulus.

**Hyperelastic advanced material**

This material class allows the user to input their own energy density function, as seen for the Neo-Hookean materials in the last section.
This can be achieved by using the extra input parameter W='f(C)' in the input file. This material needs :math:`\kappa` and :math:`\mu` as input data.
The energy density function must be a function of the form

.. math::
   W(\mu, \kappa, \mathbf{C}, \mathbf{a}) = W(\mu, \kappa, I_1, J, \mathbf{a})

with :math:`I_1=\text{trace}(\mathbf{C})`, :math:`J=\det(\mathbf{F})` and the advanced parameters :math:`\mathbf{a}` referred to as a[i] for the i-th
parameter in the function in the input file.

.. note::
   This material class uses the `num-dual <https://github.com/itt-ustutt/num-dual>`_-package for automatic differentiation. This package needs to be installed before usage:

   ``pip install num_dual``

.. autoclass:: edelweissfe.materials.hyperelasticadvanced.hyperelasticadvanced.HyperelasticAdvancedMaterial
   :members:

**Hyperelastic advanced I2 extended material**

This material class allows the user to input their own energy density function. Compared to the ``HyperelasticAdvancedMaterial``, this material allows
the second invariant :math:`I_2` and :math:`\mathbf{C}` itself to be used in the energy density function.
This can be achieved by using the extra input parameter W='f(C)' in the input file. This material needs :math:`\kappa` and :math:`\mu` as input data.
The energy density function must be a function of the form

.. math::
   W(\mu, \kappa, \mathbf{C}, \mathbf{a}) = W(\mu, \kappa, \mathbf{C}, I_1, I_2, J, \mathbf{a})

with :math:`I_1=\text{trace}(\mathbf{C})`, :math:`I_2=\frac{1}{2}\left((\text{trace}(\mathbf{C}))^2 - \text{trace}(\mathbf{C}^2)\right)`,
:math:`J=\det(\mathbf{F})` and the advanced parameters :math:`\mathbf{a}` referred to as a[i] for the i-th parameter in the function in the input file.
The right Cauchy-Green tensor itself may also be used as an input.

.. note::
   This material class uses the `autograd <https://github.com/HIPS/autograd/>`_ and `num-dual <https://github.com/itt-ustutt/num-dual>`_-package for automatic differentiation
   depending on how complicated the energy density function is. These two packages need to be installed before usage:

   ``pip install num_dual && pip install autograd``

.. autoclass:: edelweissfe.materials.hyperelasticadvanced.hyperelasticadvancedi2extended.HyperelasticAdvancedI2ExtendedMaterial
   :members:

Elastic-plastic Neo-Hookean W(I1, J) Pence-Gou [a] materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These materials use the same hyperelastic formulations as for the :ref:`Elastic Neo-Hookean W(I1, J) Pence-Gou [a] materials` with plasticity using
the spatial Hencky-strain :math:`\mathbf{e}` and the same yield and hardening functions as for the :ref:`Von Mises material` added onto them.
This material needs the following material properties in the correct order:

#. :math:`\mathbf{\mu}`          - Shear modulus.
#. :math:`\mathbf{\kappa}`       - Bulk modulus.
#. :math:`\mathbf{f_{y0}}`       - Yield stress.
#. :math:`\mathbf{H_{lin}}`      - Linear plastic hardening parameter.
#. :math:`\mathbf{\Delta f_y}`   - Multiplicator for nonlinear isotropic hardening.
#. :math:`\mathbf{\delta}`       - Exponent for nonlinear isotropic hardening.

These three material formulations solve for the residual:

.. math::
   \mathbf{R}(\bar{\mathbf{x}}) = \begin{bmatrix}\left(\mathbf{e}_{n+1}^\text{e}\right)_\text{v} - \left(\mathbf{e}_{n+1}^\text{e,trial}\right)_\text{v} +
   \sqrt{\frac{3}{2}}\Delta\bar{\kappa}\frac{\left(\tau^\text{dev}\right)_\text{v}}{||\tau^\text{dev}||} \\
   ||\tau^\text{dev}(\tau(\mathbf{e}^\text{e}_{n+1}))|| - \sqrt{\frac{2}{3}} f_\text{y}(\bar{\kappa}^\text{old}+\Delta\bar{\kappa}) \end{bmatrix}

and give back the Kirchhoff stress and its corresponding tangent modulus. The full algorithm for these materials can be found in [b], chapter 14.

.. autoclass:: edelweissfe.materials.neohookeplastic.neohookepencegouformulationaplastic.NeoHookeanWaPlasticMaterial
   :members:

.. autoclass:: edelweissfe.materials.neohookeplastic.neohookepencegouformulationbplastic.NeoHookeanWbPlasticMaterial
   :members:

.. autoclass:: edelweissfe.materials.neohookeplastic.neohookepencegouformulationcplastic.NeoHookeanWcPlasticMaterial
   :members:

Advanced hyperelastic-plastic material
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This material class allows the user to input their own energy density function.
This can be achieved by using the extra input parameter W='f(C)' in the input file. This material needs the following material properties in the correct order
The algorithm from the section :ref:`Elastic-plastic Neo-Hookean W(I1, J) Pence-Gou [a] materials` is used.
For this material, only the elastic part is differentiated by using automatic differentiation, the plastic part still uses analytical derivatives.
This material needs the following parameters as input data using keywords:

#. ``mu``       - Shear modulus.
#. ``K``        - Bulk modulus.
#. ``fy0``      - Yield stress.
#. ``HLin``     - Linear plastic hardening parameter.
#. ``dfy``      - Multiplicator for nonlinear isotropic hardening.
#. ``delta``    - Exponent for nonlinear isotropic hardening.
#. ``psi_e``    - Energy density function.
#. ``a``        - (Optional) your own parameters in ``psi_e`` (assignment with commas and quotation mark: ``a = '1.,2.,3'``).

The energy density function must be a function of the form

.. math::
   W(\mu, \kappa, I_1, J, \mathbf{a})

with :math:`\mathbf{C}` itself and :math:`I_2` being not allowed.

.. note::
   This material class uses the `num-dual <https://github.com/itt-ustutt/num-dual>`_-package for automatic differentiation. This package needs to be installed before usage:

   ``pip install num_dual``

.. autoclass:: edelweissfe.materials.hyperplasticadvanced.hyperplasticadvanced.HyperplasticAdvancedMaterial
   :members:

Implementing your own materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Relevant module for hypoelastic materials: ``edelweissfe.materials.base.basehypoelasticmaterial``

.. autoclass:: edelweissfe.materials.base.basehypoelasticmaterial.BaseHypoElasticMaterial
   :members:

Relevant module for hyperelastic materials: ``edelweissfe.materials.base.basehyperelasticmaterial``

.. autoclass:: edelweissfe.materials.base.basehyperelasticmaterial.BaseHyperElasticMaterial
   :members:


[a] Thomas J. Pence and Kun Gou, “On compressible versions of the incompressible neo-Hookean material”, Mathematics and Mechanics of Solids, 20(2): 157-182, 2015

[b] EA de Souza Neto, D Perić and DRJ Owen, ”Computational methods for plasticity - Theory and applications”, John Wiley & Sons Ltd, Engineering, Swansea University Bay Campus, Fabian Way, Swansea, SA1 8EN, 2008
