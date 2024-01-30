:py:mod:`bmlite.materials`
==========================

.. py:module:: bmlite.materials

.. autoapi-nested-parse::

   .. rubric:: Material Properties Package

   The ``materials`` package contains kinetic and transport properties for common
   battery materials, including both electrode active materials and electrolytes.



Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   bmlite.materials.Gen2Electrolyte
   bmlite.materials.GraphiteFast
   bmlite.materials.GraphiteSlow
   bmlite.materials.NMC532Fast
   bmlite.materials.NMC532Slow




.. py:class:: Gen2Electrolyte




   
   Gen 2 Electrolyte transport and material properties.

   .. py:method:: get_D(C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the lithium ion diffusivity in the electrolyte solution at
      concentration ``C_Li`` and temperature ``T``.

      :param C_Li: Lithium ion concentration in the electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **D** (*float | 1D array*) -- Lithium ion diffusivity in the electrolyte [m^2/s].


   .. py:method:: get_gamma(C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the electrolyte thermodynamic factor at concentration
      ``C_Li`` and temperature ``T``.

      :param C_Li: Lithium ion concentration in the electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **gamma** (*float | 1D array*) -- Thermodynamic factor [-].


   .. py:method:: get_kappa(C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the electrolyte conductivity at concentration ``C_Li`` and
      temperature ``T``.

      :param C_Li: Lithium ion concentration in the electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **kappa** (*float | 1D array*) -- Electrolyte conductivity [S/m].


   .. py:method:: get_t0(C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the lithium ion transference number at concentration ``C_Li``
      and temperature ``T``.

      :param C_Li: Lithium ion concentration in the electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **t0** (*float | 1D array*) -- Lithium ion transference number [-].



.. py:class:: GraphiteFast(alpha_a: float, alpha_c: float, Li_max: float)




   
   Computationally fast Graphite kinetic and transport properties.

   Differs from ``GraphiteSlow`` because the equilibrium potential is
   not piecewise here, making it less accurate, but faster to evaluate.

   :param alpha_a: Anodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_a: float
   :param alpha_c: Cathodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_c: float
   :param Li_max: Maximum lithium concentration in solid phase [kmol/m^3].
   :type Li_max: float

   .. py:method:: get_Ds(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the lithium diffusivity in the solid phase given the local
      intercalation fraction ``x`` and temperature ``T``.

      :param x: Lithium intercalation fraction [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Ds** (*float | 1D array*) -- Lithium diffusivity in the solid phase [m^2/s].


   .. py:method:: get_Eeq(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the equilibrium potential given the intercalation fraction
      ``x`` at the particle surface and temperature ``T``.

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Eeq** (*float | 1D array*) -- Equilibrium potential [V].


   .. py:method:: get_i0(x: float | numpy.ndarray, C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the exchange current density given the intercalation
      fraction ``x`` at the particle surface, the local lithium ion
      concentration ``C_Li``, and temperature ``T``. The input types for
      ``x`` and ``C_Li`` should both be the same (i.e., both float or both
      1D arrays).

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param C_Li: Lithium ion concentration in the local electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **i0** (*float | 1D array*) -- Exchange current density [A/m^2].



.. py:class:: GraphiteSlow(alpha_a: float, alpha_c: float, Li_max: float)




   
   Computationally slow Graphite kinetic and transport properties.

   Differs from ``GraphiteFast`` because the equilibrium potential is
   piecewise here, making it more accurate, but slower to evaluate.

   :param alpha_a: Anodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_a: float
   :param alpha_c: Cathodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_c: float
   :param Li_max: Maximum lithium concentration in solid phase [kmol/m^3].
   :type Li_max: float

   .. py:method:: get_Ds(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the lithium diffusivity in the solid phase given the local
      intercalation fraction ``x`` and temperature ``T``.

      :param x: Lithium intercalation fraction [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Ds** (*float | 1D array*) -- Lithium diffusivity in the solid phase [m^2/s].


   .. py:method:: get_Eeq(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the equilibrium potential given the intercalation fraction
      ``x`` at the particle surface and temperature ``T``.

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Eeq** (*float | 1D array*) -- Equilibrium potential [V].

      :raises ValueError :: x is out of bounds [x_min, x_max].


   .. py:method:: get_i0(x: float | numpy.ndarray, C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the exchange current density given the intercalation
      fraction ``x`` at the particle surface, the local lithium ion
      concentration ``C_Li``, and temperature ``T``. The input types for
      ``x`` and ``C_Li`` should both be the same (i.e., both float or both
      1D arrays).

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param C_Li: Lithium ion concentration in the local electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **i0** (*float | 1D array*) -- Exchange current density [A/m^2].



.. py:class:: NMC532Fast(alpha_a: float, alpha_c: float, Li_max: float)




   
   Computationally fast NMC532 kinetic and transport properties.

   Differs from ``NMC532Slow`` because the equilibrium potential is not
   piecewise here, making it less accurate, but faster to evaluate.

   :param alpha_a: Anodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_a: float
   :param alpha_c: Cathodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_c: float
   :param Li_max: Maximum lithium concentration in solid phase [kmol/m^3].
   :type Li_max: float

   .. py:method:: get_Ds(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the lithium diffusivity in the solid phase given the local
      intercalation fraction ``x`` and temperature ``T``.

      :param x: Lithium intercalation fraction [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Ds** (*float | 1D array*) -- Lithium diffusivity in the solid phase [m^2/s].


   .. py:method:: get_Eeq(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the equilibrium potential given the intercalation fraction
      ``x`` at the particle surface and temperature ``T``.

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Eeq** (*float | 1D array*) -- Equilibrium potential [V].


   .. py:method:: get_i0(x: float | numpy.ndarray, C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the exchange current density given the intercalation
      fraction ``x`` at the particle surface, the local lithium ion
      concentration ``C_Li``, and temperature ``T``. The input types for
      ``x`` and ``C_Li`` should both be the same (i.e., both float or both
      1D arrays).

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param C_Li: Lithium ion concentration in the local electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **i0** (*float | 1D array*) -- Exchange current density [A/m^2].



.. py:class:: NMC532Slow(alpha_a: float, alpha_c: float, Li_max: float)




   
   Computationally slow NMC532 kinetic and transport properties.

   Differs from ``NMC532Fast`` because the equilibrium potential is
   piecewise here, making it more accurate, but slower to evaluate.

   :param alpha_a: Anodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_a: float
   :param alpha_c: Cathodic symmetry factor in Butler-Volmer expression [-].
   :type alpha_c: float
   :param Li_max: Maximum lithium concentration in solid phase [kmol/m^3].
   :type Li_max: float

   .. py:method:: get_Ds(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the lithium diffusivity in the solid phase given the local
      intercalation fraction ``x`` and temperature ``T``.

      :param x: Lithium intercalation fraction [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Ds** (*float | 1D array*) -- Lithium diffusivity in the solid phase [m^2/s].


   .. py:method:: get_Eeq(x: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the equilibrium potential given the intercalation fraction
      ``x`` at the particle surface and temperature ``T``.

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **Eeq** (*float | 1D array*) -- Equilibrium potential [V].

      :raises ValueError :: x is out of bounds [x_min, x_max].


   .. py:method:: get_i0(x: float | numpy.ndarray, C_Li: float | numpy.ndarray, T: float) -> float | numpy.ndarray

      Calculate the exchange current density given the intercalation
      fraction ``x`` at the particle surface, the local lithium ion
      concentration ``C_Li``, and temperature ``T``. The input types for
      ``x`` and ``C_Li`` should both be the same (i.e., both float or both
      1D arrays).

      :param x: Lithium intercalation fraction at ``r = R_s`` [-].
      :type x: float | 1D array
      :param C_Li: Lithium ion concentration in the local electrolyte [kmol/m^3].
      :type C_Li: float | 1D array
      :param T: Battery temperature [K].
      :type T: float

      :returns: **i0** (*float | 1D array*) -- Exchange current density [A/m^2].



