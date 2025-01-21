bmlite.math
===========

.. py:module:: bmlite.math

.. autoapi-nested-parse::

   .. rubric:: Math Module

   The math module defines common mathematical operators, e.g., the gradient and
   divergence operators. Functions in this module are built to operate using
   arrays. Generally speaking, using vectorized math provides significant speed
   boosts in computational modeling by removing slower ``for`` loops and ``if``
   statements.



Functions
---------

.. autoapisummary::

   bmlite.math.div_r
   bmlite.math.div_x
   bmlite.math.grad_r
   bmlite.math.grad_x
   bmlite.math.int_r
   bmlite.math.int_x
   bmlite.math.param_combinations


Module Contents
---------------

.. py:function:: div_r(rm, rp, f, axis = -1)

   Return spherical r-divergence.

   :param rm: Independent variable values at "minus" boundaries.
   :type rm: 1D np.array
   :param rp: Independent variable values at "plus" boundaries.
   :type rp: 1D np.array
   :param f: Dependent variable evaluated at r boundaries.
   :type f: np.array
   :param axis: f dimension corresponding to r. The default is -1.
   :type axis: int, optional

   :returns: **df_dr** (*np.array*) -- The divergence ``1/r**2 * d(r**2 * f)/dr``. The shape is one fewer
             than ``f`` along the specified axis.


.. py:function:: div_x(xm, xp, f, axis = -1)

   Return Cartesian x-divergence.

   :param xm: Independent variable values at "minus" boundaries.
   :type xm: 1D np.array
   :param xp: Independent variable values at "plus" boundaries.
   :type xp: 1D np.array
   :param f: Dependent variable evaluated at x boundaries.
   :type f: np.array
   :param axis: f dimension corresponding to x. The default is -1.
   :type axis: int, optional

   :returns: **df_dx** (*np.array*) -- The divergence ``df/dx``. The shape is one fewer than ``f`` along the
             specified axis.

   .. rubric:: Notes

   This function is valid for any Cartesian direction, not just x. To get
   the 2D or 3D divergence, you have to evaluate and add the separate terms.
   For example,

   .. code-block:: python

       div_f = div_x(xm, xp, fx, 0) + div_x(ym, yp, fy, 1)
             + div_x(zm, zp, fz, 2)

   The ``fx``, ``fy``, and ``fz`` terms must be evaluated at the ``x``, ``y``,
   and ``z`` boundaries, respectively. For example, a grid with (Nx, Ny, Nz)
   volume discretizations will have an ``fx`` with shape (Nx + 1, Ny, Nz)
   because the number of boundaries is always one greater than the number of
   volumes.


.. py:function:: grad_r(r, f, axis = -1)

   Return spherical r-gradient.

   :param r: Independent variable values.
   :type r: 1D np.array
   :param f: Dependent variable values.
   :type f: np.array
   :param axis: f dimension corresponding to r. The default is -1.
   :type axis: int, optional

   :returns: **df_dr** (*np.array*) -- The gradient ``df/dr``. The shape is one fewer than ``f`` along the
             specified axis.


.. py:function:: grad_x(x, f, axis = -1)

   Return Cartesian x-gradient.

   :param x: Independent variable values.
   :type x: 1D np.array
   :param f: Dependent variable values.
   :type f: np.array
   :param axis: f dimension corresponding to x. The default is -1.
   :type axis: int, optional

   :returns: **df_dx** (*np.array*) -- The gradient ``df/dx``. The shape is one fewer than ``f`` along the
             specified axis.

   .. rubric:: Notes

   This function is valid for any Cartesian direction, not just x.


.. py:function:: int_r(rm, rp, f, axis = -1)

   Return spherical r-integral.

   :param rm: Independent variable values at "minus" boundaries.
   :type rm: 1D np.array
   :param rp: Independent variable values at "plus" boundaries.
   :type rp: 1D nb.array
   :param f: Dependent variable evaluated at r centers.
   :type f: np.array
   :param axis: f dimension corresponding to r. The default is -1.
   :type axis: int, optional

   :returns: **int_r** (*np.array*) -- The result of integration. The dimension of the result is one fewer
             than ``f`` along the specified axis.

   .. rubric:: Notes

   The result is over all spherical dimensions (r, theta, phi) assuming ``f``
   is independent of theta and phi.

   The integral is written for numerical results from finite volume solutions.
   Integration is performed over meshed control volumes, where ``f[i]`` is
   assumed uniform within a volume defined by ``rm[i] < r < rp[i]``.


.. py:function:: int_x(xm, xp, f, axis = -1)

   Return Cartesian x-integral.

   :param xm: Independent variable values at "minus" boundaries.
   :type xm: 1D np.array
   :param xp: Independent variable values at "plus" boundaries.
   :type xp: 1D np.array
   :param f: Dependent variable evaluated at x centers.
   :type f: np.array
   :param axis: f dimension corresponding to x. The default is -1.
   :type axis: int, optional

   :returns: **int_x** (*np.array*) -- The result of integration. The dimension of the result is one fewer
             than ``f`` along the specified axis.

   .. rubric:: Notes

   The integral is written for numerical results from finite volume solutions.
   Integration is performed over meshed control volumes, where ``f[i]`` is
   assumed uniform within a volume defined by ``xm[i] < x < xp[i]``.


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


