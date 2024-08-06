Materials
=========

Relevant module: ``edelweissfe.config.materiallibrary``

.. automodule:: edelweissfe.config.materiallibrary
   :members:

Materials are defined in EdelweissFE using the ``*material`` keyword.
Arguments are

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
      - (Optional) The material provider.

If the material provider is not given ``marmotmaterial`` is assumed. Materials are assigned to elements by means of :ref:`sections`.

Provider ``marmotmaterial``
---------------------------

Relevant module: `Marmot <https://github.com/MAteRialMOdelingToolbox/Marmot/>`_.

Provider ``edelweissmaterial``
------------------------------

Relevant module: ``edelweissfe.materials``

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

Implementing your own materials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Relevant module: ``edelweissfe.materials.base.basehypoelasticmaterial``

.. autoclass:: edelweissfe.materials.base.basehypoelasticmaterial.BaseHypoElasticMaterial
   :members: