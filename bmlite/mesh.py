"""
Mesh Module
-----------
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
"""

from numpy import ndarray as _ndarray


def x_ptr(domain: object, keys: list[str]) -> dict:
    """
    Add an ``x`` pointer dictionary to the domain object.

    An ``x`` pointer is a dictionary where the keys are variables and the
    values are 1D ``int`` arrays that specify the solution vector indices for
    all of the specified variable's values in the ``x`` direction. Using this
    pointer removes for loops from the model DAEs and makes it easier to
    vectorize expressions.

    Parameters
    ----------
    domain : domain object
        A domain (anode, separator, cathode, etc.) object from one of the
        model domain modules. The domain should already have ``ptr`` and
        ``Nx`` attributes. The ``ptr`` attribute is a dictionary where the
        keys are variable names and the values are the indices for the first
        occurrence of that variable. ``Nx`` is the number of ``x`` control
        volumes.

    keys : list[str]
        A list of the variable names defined in the ``x`` direction. This
        list should be a subset of the domain's existing ``ptr`` keys.

    Returns
    -------
    None.
    """

    from numpy import array

    domain.x_ptr = {}
    for k in keys:
        domain.x_ptr[k] = array([domain.ptr[k] + i * domain.ptr['x_off']
                                 for i in range(domain.Nx)], dtype=int)


def r_ptr(domain: object, keys: list[str]) -> None:
    """
    Add an ``r`` pointer dictionary to the domain object.

    An ``r`` pointer is a dictionary where the keys are variables and the
    values are 1D ``int`` arrays that specify the solution vector indices for
    all of the specified variable's values in the ``r`` direction. Using this
    pointer removes for loops from the model DAEs and makes it easier to
    vectorize expressions.

    Parameters
    ----------
    domain : domain object
        A domain (anode, separator, cathode, etc.) object from one of the
        model domain modules. The domain should already have ``ptr`` and
        ``Nr`` attributes. The ``ptr`` attribute is a dictionary where the
        keys are variable names and the values are the indices for the first
        occurrence of that variable. ``Nr`` is the number of ``r`` control
        volumes.

    keys : list[str]
        A list of the variable names defined in the ``r`` direction. This
        list should be a subset of the domain's existing ``ptr`` keys.

    Returns
    -------
    None.
    """

    from numpy import array

    domain.r_ptr = {}
    for k in keys:
        domain.r_ptr[k] = array([domain.ptr[k] + j * domain.ptr['r_off']
                                 for j in range(domain.Nr)], dtype=int)


def xr_ptr(domain: object, keys: list[str]) -> None:
    """
    Add an ``xr`` pointer dictionary to the domain object.

    An ``xr`` pointer is a dictionary where the keys are variables and the
    values are 2D ``int`` arrays that specify the solution vector indices for
    all of the specified variable's values in the ``x`` (rows) and ``r`` (cols)
    directions. Using this pointer removes for loops from the model DAEs and
    makes it easier to vectorize expressions.

    Parameters
    ----------
    domain : domain object
        A domain (anode, separator, cathode, etc.) object from one of the
        model domain modules. The domain should already have ``ptr``, ``Nx``,
        and ``Nr`` attributes. The ``ptr`` attribute is a dictionary where the
        keys are variable names and the values are the indices for the first
        occurrence of that variable. ``Nx`` and ``Nr`` are the number of ``x``
        and ``r`` control volumes, respectively.

    keys : list[str]
        A list of the variable names defined in both the ``x`` and ``r``
        directions. This list should be a subset of the domain's existing
        ``ptr`` keys.

    Returns
    -------
    None.
    """

    from numpy import array, tile

    domain.xr_ptr = {}
    for k in keys:
        first_row = array([domain.ptr[k] + j * domain.ptr['r_off']
                           for j in range(domain.Nr)], dtype=int)

        domain.xr_ptr[k] = tile(first_row, (domain.Nx, 1))

        for i in range(domain.Nx):
            domain.xr_ptr[k][i, :] += i * domain.ptr['x_off']


def uniform_mesh(Lx: float, Nx: int, x0: float = 0.) -> tuple[_ndarray]:
    """
    Determine the interface and center locations for uniformly meshed control
    volumes.

    Parameters
    ----------
    Lx : float
        The length of the domain to mesh.

    Nx : int
        Number of control volumes that span the domain's length.

    x0 : int, optional
        The reference value for the start of the domain's mesh. The default
        is 0.

    Returns
    -------
    xm : 1D array
        The "minus" interface locations. Note that ``xm[0] = x0``.

    xp : 1D array
        The "plus" interface locations. Note that ``xp[0] = x0 + Lx``.

    x : 1D array
        The locations of the control volume centers. ``xm[i] < x[i] < xp[i]``
        for all control volumes ``i``.

    Notes
    -----
    Although the inputs and outputs are labeled using ``x`` as the variable,
    this mesh is equally valid for other directions (e.g., radial).
    """

    from numpy import linspace

    xm = x0 + linspace(0., Lx * (1 - 1 / Nx), Nx)
    xp = x0 + linspace(Lx / Nx, Lx, Nx)
    x = 0.5 * (xm + xp)

    return xm, xp, x


def param_weights(xm: _ndarray, xp: _ndarray) -> tuple[_ndarray]:
    """
    Determine relative parameter weights between adjacent control volumes.

    In finite volume method solutions, fluxes are calculated at the boundaries
    between control volumes. Therefore, locally resolved parameters must be
    weighted between the adjacent volumes participating in the calculation.
    This function calculates those weights.

    Parameters
    ----------
    xm : 1D array
        Locations of the "minus" interfaces in the mesh.

    xp : 1D array
        Locations of the "plus" interfaces in the mesh.

    Returns
    -------
    wt_m : 1D array
        Parameter weights for each of the "minus" half volumes.

    wt_p : 1D array
        Parameter weights for each of the "plus" half volumes.
    """

    x = 0.5 * (xm + xp)

    wt_m = 0.5 * (xp[:-1] - xm[:-1]) / (x[1:] - x[:-1])
    wt_p = 0.5 * (xp[1:] - xm[1:]) / (x[1:] - x[:-1])

    return wt_m, wt_p
