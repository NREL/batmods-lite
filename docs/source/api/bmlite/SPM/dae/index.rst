bmlite.SPM.dae
==============

.. py:module:: bmlite.SPM.dae

.. autoapi-nested-parse::

   .. rubric:: DAE Module

   This module includes the system of differential algebraic equations (DAE) for
   the SPM model. In addition, the ``bandwidth`` function is included in this
   module, which helps determine the lower and upper bandwidths of ``residuals``
   so the ``'band'`` linear solver option can be used in the ``IDASolver`` class.



Functions
---------

.. autoapisummary::

   bmlite.SPM.dae.residuals


Module Contents
---------------

.. py:function:: residuals(t, sv, svdot, res, inputs)

   The DAE residuals ``res = M*y' - f(t, y)`` for the SPM model.

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
   :type inputs: (sim : SPM Simulation object, exp : experiment dict)

   :returns: **outputs** (*tuple[np.ndarray]*) -- If the experimental step `mode` is set to `post`, then the following
             post-processed variables will be returned in a tuple. Otherwise,
             returns None.

             ========= =================================================
             Variable  Description [units] (*type*)
             ========= =================================================
             sdot_an   anode Li+ production [kmol/m^3/s] (*1D array*)
             sdot_ca   cathode Li+ production [kmol/m^3/s] (*1D array*)
             ========= =================================================


