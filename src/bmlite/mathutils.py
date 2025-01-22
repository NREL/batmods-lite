"""
Math Module
-----------
The math module defines common mathematical operators, e.g., the gradient and
divergence operators. Functions in this module are built to operate using
arrays. Generally speaking, using vectorized math provides significant speed
boosts in computational modeling by removing slower ``for`` loops and ``if``
statements.

"""

import numpy as np


def grad_x(x: np.ndarray, f: np.ndarray, axis: int = -1) -> np.ndarray:
    """
    Return Cartesian x-gradient.

    Parameters
    ----------
    x : 1D np.array
        Independent variable values.
    f : np.array
        Dependent variable values.
    axis : int, optional
        f dimension corresponding to x. The default is -1.

    Returns
    -------
    df_dx : np.array
        The gradient ``df/dx``. The shape is one fewer than ``f`` along the
        specified axis.

    Notes
    -----
    This function is valid for any Cartesian direction, not just x.

    """

    new_axis = [1] * f.ndim
    new_axis[axis] = -1

    dx = x[1:] - x[:-1]

    df_dx = np.diff(f, axis=axis) / dx.reshape(new_axis)

    return df_dx


def grad_r(r: np.ndarray, f: np.ndarray, axis: int = -1) -> np.ndarray:
    """
    Return spherical r-gradient.

    Parameters
    ----------
    r : 1D np.array
        Independent variable values.
    f : np.array
        Dependent variable values.
    axis : int, optional
        f dimension corresponding to r. The default is -1.

    Returns
    -------
    df_dr : np.array
        The gradient ``df/dr``. The shape is one fewer than ``f`` along the
        specified axis.

    """

    new_axis = [1] * f.ndim
    new_axis[axis] = -1

    dr = r[1:] - r[:-1]

    df_dr = np.diff(f, axis=axis) / dr.reshape(new_axis)

    return df_dr


def div_x(xm: np.ndarray, xp: np.ndarray, f: np.ndarray,
          axis: int = -1) -> np.ndarray:
    """
    Return Cartesian x-divergence.

    Parameters
    ----------
    xm : 1D np.array
        Independent variable values at "minus" boundaries.
    xp : 1D np.array
        Independent variable values at "plus" boundaries.
    f : np.array
        Dependent variable evaluated at x boundaries.
    axis : int, optional
        f dimension corresponding to x. The default is -1.

    Returns
    -------
    df_dx : np.array
        The divergence ``df/dx``. The shape is one fewer than ``f`` along the
        specified axis.

    Notes
    -----
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

    """

    new_axis = [1] * f.ndim
    new_axis[axis] = -1

    dx = xp - xm

    df_dx = np.diff(f, axis=axis) / dx.reshape(new_axis)

    return df_dx


def div_r(rm: np.ndarray, rp: np.ndarray, f: np.ndarray,
          axis: int = -1) -> np.ndarray:
    """
    Return spherical r-divergence.

    Parameters
    ----------
    rm : 1D np.array
        Independent variable values at "minus" boundaries.
    rp : 1D np.array
        Independent variable values at "plus" boundaries.
    f : np.array
        Dependent variable evaluated at r boundaries.
    axis : int, optional
        f dimension corresponding to r. The default is -1.

    Returns
    -------
    df_dr : np.array
        The divergence ``1/r**2 * d(r**2 * f)/dr``. The shape is one fewer
        than ``f`` along the specified axis.

    """

    new_axis = [1] * f.ndim
    new_axis[axis] = -1

    rm = rm.reshape(new_axis)
    rp = rp.reshape(new_axis)

    r = 0.5 * (rm + rp)
    dr = rp - rm

    fm = np.delete(f, -1, axis=axis)
    fp = np.delete(f, 0, axis=axis)

    df_dr = 1. / r**2 * (rp**2 * fp - rm**2 * fm) / dr

    return df_dr


def int_x(xm: np.ndarray, xp: np.ndarray, f: np.ndarray,
          axis: int = -1) -> np.ndarray:
    """
    Return Cartesian x-integral.

    Parameters
    ----------
    xm : 1D np.array
        Independent variable values at "minus" boundaries.
    xp : 1D np.array
        Independent variable values at "plus" boundaries.
    f : np.array
        Dependent variable evaluated at x centers.
    axis : int, optional
        f dimension corresponding to x. The default is -1.

    Returns
    -------
    int_x : np.array
        The result of integration. The dimension of the result is one fewer
        than ``f`` along the specified axis.

    Notes
    -----
    The integral is written for numerical results from finite volume solutions.
    Integration is performed over meshed control volumes, where ``f[i]`` is
    assumed uniform within a volume defined by ``xm[i] < x < xp[i]``.

    """

    new_axis = [1] * f.ndim
    new_axis[axis] = -1

    xm = xm.reshape(new_axis)
    xp = xp.reshape(new_axis)

    int_x = np.sum(f * (xp - xm), axis=axis)

    return int_x


def int_r(rm: np.ndarray, rp: np.ndarray, f: np.ndarray,
          axis: int = -1) -> np.ndarray:
    """
    Return spherical r-integral.

    Parameters
    ----------
    rm : 1D np.array
        Independent variable values at "minus" boundaries.
    rp : 1D nb.array
        Independent variable values at "plus" boundaries.
    f : np.array
        Dependent variable evaluated at r centers.
    axis : int, optional
        f dimension corresponding to r. The default is -1.

    Returns
    -------
    int_r : np.array
        The result of integration. The dimension of the result is one fewer
        than ``f`` along the specified axis.

    Notes
    -----
    The result is over all spherical dimensions (r, theta, phi) assuming ``f``
    is independent of theta and phi.

    The integral is written for numerical results from finite volume solutions.
    Integration is performed over meshed control volumes, where ``f[i]`` is
    assumed uniform within a volume defined by ``rm[i] < r < rp[i]``.

    """

    new_axis = [1] * f.ndim
    new_axis[axis] = -1

    rm = rm.reshape(new_axis)
    rp = rp.reshape(new_axis)

    r = 0.5 * (rm + rp)

    int_r = np.sum(4. * np.pi * r**2 * f * (rp - rm), axis=axis)

    return int_r


def param_combinations(params: list[str],
                       values: list[np.ndarray]) -> list[dict]:
    """
    Generate all possible combinations for a set of parameters given their
    possible values.

    Parameters
    ----------
        params : list[str]
            List of parameter names, including the domain, e.g. ``an.i0_deg``.
        values : list[1D array]
            List of possible values for each parameter. The array in each
            index ``i`` should correspond to the variable given in ``params``
            with the same index.

    Returns
    -------
        combinations : list[dict]
            List of dictionaries representing all possible combinations of
            parameter values.

    """

    from itertools import product

    combinations = []
    for combination in product(*values):
        combinations.append({k: v for k, v in zip(params, combination)})

    return combinations
