:py:mod:`bmlite.mesh`
=====================

.. py:module:: bmlite.mesh

.. autoapi-nested-parse::

   .. rubric:: Mesh Module

   This module contains functions to build pointers and meshes. Pointers are
   the integer indices for state variables. State variables in a model may be
   defined in different dimensions. For example, ``phi_ed`` is 0D in the single
   particle model and ``Li_ed`` is 1D in the ``r`` direction. Therefore, within
   any given model, it may be useful to have a combination of ``ptr`` for 0D,
   ``x_ptr`` for 1D in ``x``, ``xr_ptr`` for 2D in ``x`` and ``r``, etc.

   This module also contains functions to build meshes. In finite volume methods,
   the mesh specifies the locations of the control volume centers as well as the
   interfaces between two control volumes. It makes sense to store both pointers
   and mesh functions in the same module because the pointers must be self
   consistent for any given model.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bmlite.mesh.param_weights
   bmlite.mesh.r_ptr
   bmlite.mesh.uniform_mesh
   bmlite.mesh.x_ptr
   bmlite.mesh.xr_ptr



.. py:function:: param_weights(xm, xp)

   Determine relative parameter weights between adjacent control volumes.

   In finite volume method solutions, fluxes are calculated at the boundaries
   between control volumes. Therefore, locally resolved parameters must be
   weighted between the adjacent volumes participating in the calculation.
   This function calculates those weights.

   :param xm: Locations of the "minus" interfaces in the mesh.
   :type xm: 1D array
   :param xp: Locations of the "plus" interfaces in the mesh.
   :type xp: 1D array

   :returns: * **wt_m** (*1D array*) -- Parameter weights for each of the "minus" half volumes.
             * **wt_p** (*1D array*) -- Parameter weights for each of the "plus" half volumes.


.. py:function:: r_ptr(domain, keys)

   Add an ``r`` pointer dictionary to the domain object.

   An ``r`` pointer is a dictionary where the keys are variables and the
   values are 1D ``int`` arrays that specify the solution vector indices for
   all of the specified variable's values in the ``r`` direction. Using this
   pointer removes for loops from the model DAEs and makes it easier to
   vectorize expressions.

   :param domain: A domain (anode, separator, cathode, etc.) object from one of the
                  model domain modules. The domain should already have ``ptr`` and
                  ``Nr`` attributes. The ``ptr`` attribute is a dictionary where the
                  keys are variable names and the values are the indices for the first
                  occurrence of that variable. ``Nr`` is the number of ``r`` control
                  volumes.
   :type domain: domain object
   :param keys: A list of the variable names defined in the ``r`` direction. This
                list should be a subset of the domain's existing ``ptr`` keys.
   :type keys: list[str]

   :returns: *None.*


.. py:function:: uniform_mesh(Lx, Nx, x0 = 0.0)

   Determine the interface and center locations for uniformly meshed control
   volumes.

   :param Lx: The length of the domain to mesh.
   :type Lx: float
   :param Nx: Number of control volumes that span the domain's length.
   :type Nx: int
   :param x0: The reference value for the start of the domain's mesh. The default
              is 0.
   :type x0: int, optional

   :returns: * **xm** (*1D array*) -- The "minus" interface locations. Note that ``xm[0] = x0``.
             * **xp** (*1D array*) -- The "plus" interface locations. Note that ``xp[0] = x0 + Lx``.
             * **x** (*1D array*) -- The locations of the control volume centers. ``xm[i] < x[i] < xp[i]``
               for all control volumes ``i``.

   .. rubric:: Notes

   Although the inputs and outputs are labeled using ``x`` as the variable,
   this mesh is equally valid for other directions (e.g., radial).


.. py:function:: x_ptr(domain, keys)

   Add an ``x`` pointer dictionary to the domain object.

   An ``x`` pointer is a dictionary where the keys are variables and the
   values are 1D ``int`` arrays that specify the solution vector indices for
   all of the specified variable's values in the ``x`` direction. Using this
   pointer removes for loops from the model DAEs and makes it easier to
   vectorize expressions.

   :param domain: A domain (anode, separator, cathode, etc.) object from one of the
                  model domain modules. The domain should already have ``ptr`` and
                  ``Nx`` attributes. The ``ptr`` attribute is a dictionary where the
                  keys are variable names and the values are the indices for the first
                  occurrence of that variable. ``Nx`` is the number of ``x`` control
                  volumes.
   :type domain: domain object
   :param keys: A list of the variable names defined in the ``x`` direction. This
                list should be a subset of the domain's existing ``ptr`` keys.
   :type keys: list[str]

   :returns: *None.*


.. py:function:: xr_ptr(domain, keys)

   Add an ``xr`` pointer dictionary to the domain object.

   An ``xr`` pointer is a dictionary where the keys are variables and the
   values are 2D ``int`` arrays that specify the solution vector indices for
   all of the specified variable's values in the ``x`` (rows) and ``r`` (cols)
   directions. Using this pointer removes for loops from the model DAEs and
   makes it easier to vectorize expressions.

   :param domain: A domain (anode, separator, cathode, etc.) object from one of the
                  model domain modules. The domain should already have ``ptr``, ``Nx``,
                  and ``Nr`` attributes. The ``ptr`` attribute is a dictionary where the
                  keys are variable names and the values are the indices for the first
                  occurrence of that variable. ``Nx`` and ``Nr`` are the number of ``x``
                  and ``r`` control volumes, respectively.
   :type domain: domain object
   :param keys: A list of the variable names defined in both the ``x`` and ``r``
                directions. This list should be a subset of the domain's existing
                ``ptr`` keys.
   :type keys: list[str]

   :returns: *None.*


