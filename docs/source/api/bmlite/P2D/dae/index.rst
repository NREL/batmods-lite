bmlite.P2D.dae
==============

.. py:module:: bmlite.P2D.dae

.. autoapi-nested-parse::

   .. rubric:: DAE Module

   This module includes the system of differential algebraic equations (DAE) for
   the P2D model. In addition, the ``bandwidth`` function is included in this
   module, which helps determine the lower and upper bandwidths of ``residuals``
   so the ``'band'`` linear solver option can be used in the ``IDASolver`` class.



Functions
---------

.. autoapisummary::

   bmlite.P2D.dae.residuals


Module Contents
---------------

.. py:function:: residuals(t, sv, svdot, res, inputs)

   The DAE residuals ``res = M*y' - f(t, y)`` for the P2D model.

   :param t: Value of time [s].
   :type t: float
   :param sv: Solution/state variables at time ``t``.
   :type sv: 1D array
   :param svdot: Solution/state variable time derivatives at time ``t``.
   :type svdot: 1D array
   :param res: An array the same size as ``sv`` and ``svdot``. The values are filled
               in with ``res = M*y' - f(t, y)`` inside this function.
   :type res: 1D array
   :param inputs: The simulation object and experimental details dictionary inputs that
                  describe the specific battery and experiment to simulate.
   :type inputs: (sim : P2D Simulation object, exp : experiment dict)

   :returns: **outputs** (*tuple[np.ndarray]*) -- If the experimental step `mode` is set to `post`, then the following
             post-processed variables will be returned in a tuple. Otherwise,
             returns None.

             ========== ======================================================
             Variable   Description [units] (*type*)
             ========== ======================================================
             div_i_an   divergence in anode volume [A/m3] (*1D array*)
             div_i_sep  divergence in separator volume [A/m3] (*1D array*)
             div_i_ca   divergence in cathode volume [A/m3] (*1D array*)
             sdot_an    Li+ production at each x_a [kmol/m^3/s] (*1D array*)
             sdot_ca    Li+ production at each x_c [kmol/m^3/s] (*1D array*)
             sum_ip     i_ed + i_el at "plus" interfaces [A/m^2] (*1D array*)
             i_el_x     i_el at each x interface [A/m^2] (*1D array*)
             ========== ======================================================


