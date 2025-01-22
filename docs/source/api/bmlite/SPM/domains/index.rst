bmlite.SPM.domains
==================

.. py:module:: bmlite.SPM.domains

.. autoapi-nested-parse::

   .. rubric:: Domains Module

   Contains classes to construct the battery for SPM simulations. Each class reads
   in keyword arguments that define parameters relevant to its specific domain.
   For example, the area and temperature are ``Battery`` level parameters because
   they are the same everywhere, but the discretization ``Nr`` may be different
   for the anode and cathode domains.



Classes
-------

.. autoapisummary::

   bmlite.SPM.domains.Battery
   bmlite.SPM.domains.Electrode
   bmlite.SPM.domains.Electrolyte


Module Contents
---------------

.. py:class:: Battery(**kwargs)

   
   A class for battery-level attributes.

   :param \*\*kwargs: Keyword arguments to set the battery-level attributes. The required
                      keys and descriptions are given below:

                      ====== =====================================================
                      Key    Description [units] (type)
                      ====== =====================================================
                      cap    nominal battery capacity [A*h] (*float*)
                      temp   temperature [K] (*float*)
                      area   area normal to the current collectors [m^2] (*float*)
                      ====== =====================================================
   :type \*\*kwargs: dict, required


   .. py:method:: update()

      Updates any secondary/dependent parameters. At present, this method
      does not do anything for the ``Battery`` class.

      :returns: *None.*



.. py:class:: Electrode(name, **kwargs)

   
   A class for the electrode-specific attributes and methods.

   This class is used to build both the anode and cathode for P2D
   simulations.

   :param name: Name of the electrode, must be either 'anode' or 'cathode'.
   :type name: str
   :param \*\*kwargs: Keyword arguments to set the electrode attributes. The required
                      keys and descriptions are given below:

                      ========= =========================================================
                      Key       Description [units] (type)
                      ========= =========================================================
                      Nr        number of ``r`` discretizations [-] (*int*)
                      thick     electrode thickness [m] (*float*)
                      R_s       represenatative particle radius [m] (*float*)
                      eps_el    electrolyte volume fraction [-] (*float*)
                      eps_CBD   carbon binder domain volume fraction [-] (*float*)
                      eps_void  Void volume fraction [-] (*float*)
                      alpha_a   Butler-Volmer anodic symmetry factor [-] (*float*)
                      alpha_c   Butler-Volmer cathodic symmetry factor [-] (*float*)
                      Li_max    max solid-phase lithium concentration [kmol/m^3] (*float*)
                      x_0       initial solid-phase intercalation fraction [-] (*float*)
                      i0_deg    ``i0`` degradation factor [-] (*float*)
                      Ds_deg    ``Ds`` degradation factor [-] (*float*)
                      material  class name from ``bmlite.materials`` [-] (*str*)
                      ========= =========================================================
   :type \*\*kwargs: dict, required


   .. py:method:: get_Ds(x, T)

      Calculate the lithium diffusivity in the solid phase given the local
      intercalation fraction ``x`` and temperature ``T``.

      :param x: Lithium intercalation fraction [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Ds** (*float | 1D array*) -- Lithium diffusivity in the solid phase [m^2/s].



   .. py:method:: get_Eeq(x)

      Calculate the equilibrium potential given the surface intercalation
      fraction ``x`` at the particle surface and temperature ``T``.

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float

      :returns: **Eeq** (*float*) -- Equilibrium potential [V].



   .. py:method:: get_i0(x, C_Li, T)

      Calculate the exchange current density given the surface intercalation
      fraction ``x`` at the particle surface, the local lithium ion
      concentration ``C_Li``, and temperature ``T``.

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float
      :param C_Li: Lithium ion concentration in the local electrolyte [kmol/m^3].
      :type C_Li: float
      :param T: Battery temperature [K].
      :type T: float

      :returns: **i0** (*float*) -- Exchange current density [A/m^2].



   .. py:method:: make_mesh(pshift = 0)

      Determines/sets the ``r`` locations for all of the "minus" interfaces
      ``rm``, "plus" interfaces ``rp``, and control volume centers ``r``
      based on the representative particle radius and ``Nr`` discretization.
      At present, only a uniform mesh is supported.

      :param pshift:
      :type pshift: int, optional

      :returns: *None.*

      .. seealso:: :obj:`batmods.mesh.r_ptr`, :obj:`batmods.mesh.uniform_mesh`



   .. py:method:: update()

      Updates any secondary/dependent parameters. For the ``Electrode``
      class, this initializes the material class, and sets the following:

      * Void-phase volume fraction [-]:
          ``eps_void = 1 - eps_s - eps_el``
      * Activate material volume fraction [-]:
          ``eps_AM = eps_s - eps_CBD``
      * Specific particle surface area [m^2/m^3]:
          ``A_s = 3 * eps_AM / R_s``

      :returns: *None.*



.. py:class:: Electrolyte(**kwargs)

   
   A class for the electrolyte attributes and methods.

   :param \*\*kwargs: Keyword arguments to set the electrolyte attributes. The required
                      keys and descriptions are given below:

                      ======== ==============================================
                      Key      Description [units] (type)
                      ======== ==============================================
                      li_0     initial Li+ concentration [kmol/m^3] (*float*)
                      ======== ==============================================
   :type \*\*kwargs: dict, required


   .. py:method:: update()

      Updates any secondary/dependent parameters. At present, this method
      does not do anything for the ``Electrolyte`` class.

      :returns: *None.*



