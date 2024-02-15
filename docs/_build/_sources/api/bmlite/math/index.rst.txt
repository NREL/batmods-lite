:py:mod:`bmlite.math`
=====================

.. py:module:: bmlite.math

.. autoapi-nested-parse::

   .. rubric:: Math Module

   The math module defines common mathematical operators, e.g., the gradient and
   divergence operators. Functions in this module are built to operate using
   arrays. Generally speaking, using vectorized math provides significant speed
   boosts in computational modeling by removing slower ``for`` loops and ``if``
   statements.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.math.div_r
   bmlite.math.div_x
   bmlite.math.grad_r
   bmlite.math.grad_x
   bmlite.math.int_r
   bmlite.math.int_x
   bmlite.math.param_combinations



.. py:function:: div_r(rm, rp, f, axis = 0)

   The divergence of vector field ``f`` in the spherical ``r`` direction.

   If ``f`` is a multi-dimensional array (i.e., greater than 1D), then
   ``axis`` should specify the dimension that is consistent with the
   specified ``r`` array.

   :param rm: Coordinate values for the control volumes' "minus" interfaces.
   :type rm: 1D array
   :param rp: Coordinate values for the control volumes' "plus" interfaces.
   :type rp: 1D array
   :param f: Dependent variable values that correspond to the control volumes'
             interface coordinates along ``axis``.
   :type f: ND array
   :param axis: The axis index of ``f`` that corresponds to the ``rm`` and ``rp``
                interface coordinate values. The default is 0.
   :type axis: int, optional

   :returns: **df_dr** (*ND array*) -- The divergence ``1/r**2 * d(r**2 * f)/dr``. Note that the shape of
             ``df_dr`` will be reduced by one, compared to ``f``, along ``axis``.


.. py:function:: div_x(xm, xp, f, axis = 0)

   The divergence of vector field ``f`` in the ``x`` direction.

   Although ``x`` is used here as the variable name, this function is equally
   valid for any other Cartesian coordinate. If ``f`` is a multi-dimensional
   array (i.e., greater than 1D), then ``axis`` should specify the dimension
   that is consistent with the specified ``x`` array.

   To get the total divergence for an N-dimensional array, you will have to
   add together the divergence components. For example,

   .. code-block:: python

       div_f = div_x(xm, xp, fx, 0) + div_x(ym, yp, fy, 1)
             + div_x(zm, zp, fz, 2)

   where ``fx`` is a vector field evaluated at all ``x`` interfaces that
   intersect the ``y`` and ``z`` control volume centers, ``fy`` is a vector
   field evaluated at all ``y`` interfaces that intersect the ``x`` and ``z``
   control volume centers, and ``fz`` is a vector field evaluated at all
   ``z`` interfaces that intersect the ``x`` and ``y`` control volume centers.
   Therefore, if ``(x, y, z)`` has shape ``(n, m, p)`` then ``fx``, ``fy``,
   and ``fz`` have the respective shapes ``(n+1, m, p)``, ``(n+1, m+1, p)``,
   and ``(n+1, m, p+1)``

   :param xm: Coordinate values for the control volumes' "minus" interfaces.
   :type xm: 1D array
   :param xp: Coordinate values for the control voluemes' "plus" interfaces.
   :type xp: 1D array
   :param f: Dependent variable values that correspond to the control volumes'
             interface coordinates along ``axis``.
   :type f: ND array
   :param axis: The axis index of ``f`` that corresponds to the ``xm`` and ``xp``
                interface coordinate values. The default is 0.
   :type axis: int, optional

   :returns: **df_dx** (*ND array*) -- The divergence ``df/dx``. Note that the shape of ``df_dx`` will be
             reduced by one, compared to ``f``, along ``axis``.


.. py:function:: grad_r(r, f, axis = 0)

   The gradient of ``f`` in the spherical ``r`` direction.

   If ``f`` is a multi-dimensional array (i.e., greater than 1D), then
   ``axis`` should specify the dimension that is consistent with the
   specified ``r`` array.

   :param r: Coordinate values for the independent variable.
   :type r: 1D array
   :param f: Dependent variable values that correspond to the ``r`` coordinate
             locations along ``axis``.
   :type f: ND array
   :param axis: The axis index of ``f`` that corresponds to the ``r`` coordinate
                values. The default is 0.
   :type axis: int, optional

   :returns: **df_dr** (*ND array*) -- The gradient ``df/dr``. Note that the shape of ``df_dr`` will be
             reduced by one, compared to ``f``, along ``axis``.


.. py:function:: grad_x(x, f, axis = 0)

   The gradient of ``f`` in the ``x`` direction.

   Although ``x`` is used here as the variable name, this function is equally
   valid for any other Cartesian coordinate. If ``f`` is a multi-dimensional
   array (i.e., greater than 1D), then ``axis`` should specify the dimension
   that is consistent with the specified ``x`` array.

   :param x: Coordinate values for the independent variable.
   :type x: 1D array
   :param f: Dependent variable values that correspond to the ``x`` coordinate
             locations along ``axis``.
   :type f: ND array
   :param axis: The axis index of ``f`` that corresponds to the ``x`` coordinate
                values. The default is 0.
   :type axis: int, optional

   :returns: **df_dx** (*ND array*) -- The gradient ``df/dx``. Note that the shape of ``df_dx`` will be
             reduced by one, compared to ``f``, along ``axis``.


.. py:function:: int_r(rm, rp, f, axis = 0)

   Integral of ``f`` with respect to the spherical ``r`` assuming ``f``
   is uniform within each control volume ``i`` defined by the bounds
   ``rm[i] <= r[i] <= rp[i]``.

   :param rm: Coordinate values for the control volumes' "minus" interfaces.
   :type rm: 1D array
   :param rp: Coordinate values for the control voluemes' "plus" interfaces.
   :type rp: 1D array
   :param f: Dependent variable values that correspond to the control volumes'
             center coordinates along ``axis``.
   :type f: ND array
   :param axis: The axis index of ``f`` that corresponds to the ``x`` coordinate. The
                default is 0.
   :type axis: int, optional

   :returns: **int_r** (*float | ND array*) -- The result of integration. If ``f`` is 1D, the returned value will be
             a scalar. Otherwise, an array with one fewer dimensions than ``f`` will
             be returned.


.. py:function:: int_x(xm, xp, f, axis = 0)

   Integral of ``f`` with respect to ``x`` assuming ``f`` is uniform within
   each control volume ``i`` defined by the bounds ``xm[i] <= x[i] <= xp[i]``.

   :param xm: Coordinate values for the control volumes' "minus" interfaces.
   :type xm: 1D array
   :param xp: Coordinate values for the control voluemes' "plus" interfaces.
   :type xp: 1D array
   :param f: Dependent variable values that correspond to the control volumes'
             center coordinates along ``axis``.
   :type f: ND array
   :param axis: The axis index of ``f`` that corresponds to the ``x`` coordinate. The
                default is 0.
   :type axis: int, optional

   :returns: **int_x** (*float | ND array*) -- The result of integration. If ``f`` is 1D, the returned value will be
             a scalar. Otherwise, an array with one fewer dimensions than ``f`` will
             be returned.


.. py:function:: param_combinations(params, values)

   Generate all possible combinations for a set of parameters given their
   possible values.

   :param params: List of parameter names, including the domain, e.g. ``an.i0_deg``.
   :type params: list[str]
   :param values: List of possible values for each parameter. The array in each
                  index ``i`` should correspond to the variable given in ``params``
                  with the same index.
   :type values: list[1D array]

   :returns: **combinations** (*list[dict]*) -- List of dictionaries representing all possible combinations of
             parameter values.


