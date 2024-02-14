"""
Math Module
-----------
The math module defines common mathematical operators, e.g., the gradient and
divergence operators. Functions in this module are built to operate using
arrays. Generally speaking, using vectorized math provides significant speed
boosts in computational modeling by removing slower ``for`` loops and ``if``
statements.
"""

from numpy import ndarray as _ndarray


def grad_x(x: _ndarray, f: _ndarray, axis: int = 0) -> _ndarray:
    """
    The gradient of ``f`` in the ``x`` direction.

    Although ``x`` is used here as the variable name, this function is equally
    valid for any other Cartesian coordinate. If ``f`` is a multi-dimensional
    array (i.e., greater than 1D), then ``axis`` should specify the dimension
    that is consistent with the specified ``x`` array.

    Parameters
    ----------
    x : 1D array
        Coordinate values for the independent variable.

    f : ND array
        Dependent variable values that correspond to the ``x`` coordinate
        locations along ``axis``.

    axis : int, optional
        The axis index of ``f`` that corresponds to the ``x`` coordinate
        values. The default is 0.

    Returns
    -------
    df_dx : ND array
        The gradient ``df/dx``. Note that the shape of ``df_dx`` will be
        reduced by one, compared to ``f``, along ``axis``.
    """

    from numpy import expand_dims

    size = f.shape[axis]

    dx = x[1:] - x[:-1]

    new_axes = [*range(f.ndim)]
    new_axes.remove(axis)
    dx = expand_dims(dx, axis=new_axes)

    f1 = f.take(range(0, size - 1), axis=axis)
    f2 = f.take(range(1, size), axis=axis)

    df_dx = (f2 - f1) / dx

    return df_dx


def grad_r(r: _ndarray, f: _ndarray, axis: int = 0) -> _ndarray:
    """
    The gradient of ``f`` in the spherical ``r`` direction.

    If ``f`` is a multi-dimensional array (i.e., greater than 1D), then
    ``axis`` should specify the dimension that is consistent with the
    specified ``r`` array.

    Parameters
    ----------
    r : 1D array
        Coordinate values for the independent variable.

    f : ND array
        Dependent variable values that correspond to the ``r`` coordinate
        locations along ``axis``.

    axis : int, optional
        The axis index of ``f`` that corresponds to the ``r`` coordinate
        values. The default is 0.

    Returns
    -------
    df_dr : ND array
        The gradient ``df/dr``. Note that the shape of ``df_dr`` will be
        reduced by one, compared to ``f``, along ``axis``.
    """

    from numpy import expand_dims

    size = f.shape[axis]

    dr = r[1:] - r[:-1]

    new_axes = [*range(f.ndim)]
    new_axes.remove(axis)
    dr = expand_dims(dr, axis=new_axes)

    f1 = f.take(range(0, size - 1), axis=axis)
    f2 = f.take(range(1, size), axis=axis)

    df_dr = (f2 - f1) / dr

    return df_dr


def div_x(xm: _ndarray, xp: _ndarray, f: _ndarray,
          axis: int = 0) -> _ndarray:
    """
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

    Parameters
    ----------
    xm : 1D array
        Coordinate values for the control volumes' "minus" interfaces.

    xp : 1D array
        Coordinate values for the control voluemes' "plus" interfaces.

    f : ND array
        Dependent variable values that correspond to the control volumes'
        interface coordinates along ``axis``.

    axis : int, optional
        The axis index of ``f`` that corresponds to the ``xm`` and ``xp``
        interface coordinate values. The default is 0.

    Returns
    -------
    df_dx : ND array
        The divergence ``df/dx``. Note that the shape of ``df_dx`` will be
        reduced by one, compared to ``f``, along ``axis``.
    """

    from numpy import expand_dims

    size = f.shape[axis]

    dx = xp - xm

    new_axes = [*range(f.ndim)]
    new_axes.remove(axis)
    dx = expand_dims(dx, axis=new_axes)

    f1 = f.take(range(0, size - 1), axis=axis)
    f2 = f.take(range(1, size), axis=axis)

    df_dx = (f2 - f1) / dx

    return df_dx


def div_r(rm: _ndarray, rp: _ndarray, f: _ndarray,
          axis: int = 0) -> _ndarray:
    """
    The divergence of vector field ``f`` in the spherical ``r`` direction.

    If ``f`` is a multi-dimensional array (i.e., greater than 1D), then
    ``axis`` should specify the dimension that is consistent with the
    specified ``r`` array.

    Parameters
    ----------
    rm : 1D array
        Coordinate values for the control volumes' "minus" interfaces.

    rp : 1D array
        Coordinate values for the control volumes' "plus" interfaces.

    f : ND array
        Dependent variable values that correspond to the control volumes'
        interface coordinates along ``axis``.

    axis : int, optional
        The axis index of ``f`` that corresponds to the ``rm`` and ``rp``
        interface coordinate values. The default is 0.

    Returns
    -------
    df_dr : ND array
        The divergence ``1/r**2 * d(r**2 * f)/dr``. Note that the shape of
        ``df_dr`` will be reduced by one, compared to ``f``, along ``axis``.
    """

    from numpy import expand_dims

    size = f.shape[axis]

    r = 0.5 * (rm + rp)
    dr = rp - rm

    new_axes = [*range(f.ndim)]
    new_axes.remove(axis)

    r = expand_dims(r, axis=new_axes)
    dr = expand_dims(dr, axis=new_axes)
    rm = expand_dims(rm, axis=new_axes)
    rp = expand_dims(rp, axis=new_axes)

    f1 = f.take(range(0, size - 1), axis=axis)
    f2 = f.take(range(1, size), axis=axis)

    df_dr = 1 / r**2 * (rp**2 * f2 - rm**2 * f1) / dr

    return df_dr


def int_x(xm: _ndarray, xp: _ndarray, f: _ndarray,
          axis: int = 0) -> float | _ndarray:
    """
    Integral of ``f`` with respect to ``x`` assuming ``f`` is uniform within
    each control volume ``i`` defined by the bounds ``xm[i] <= x[i] <= xp[i]``.

    Parameters
    ----------
    xm : 1D array
        Coordinate values for the control volumes' "minus" interfaces.

    xp : 1D array
        Coordinate values for the control voluemes' "plus" interfaces.

    f : ND array
        Dependent variable values that correspond to the control volumes'
        center coordinates along ``axis``.

    axis : int, optional
        The axis index of ``f`` that corresponds to the ``x`` coordinate. The
        default is 0.

    Returns
    -------
    int_x : float | ND array
        The result of integration. If ``f`` is 1D, the returned value will be
        a scalar. Otherwise, an array with one fewer dimensions than ``f`` will
        be returned.
    """

    import numpy as np

    new_axes = [*range(f.ndim)]
    new_axes.remove(axis)

    xm = np.expand_dims(xm, axis=new_axes)
    xp = np.expand_dims(xp, axis=new_axes)

    int_x = np.sum(f * (xp - xm), axis=axis)

    return int_x


def int_r(rm: _ndarray, rp: _ndarray, f: _ndarray,
          axis: int = 0) -> float | _ndarray:
    """
    Integral of ``f`` with respect to the spherical ``r`` assuming ``f``
    is uniform within each control volume ``i`` defined by the bounds
    ``rm[i] <= r[i] <= rp[i]``.

    Parameters
    ----------
    rm : 1D array
        Coordinate values for the control volumes' "minus" interfaces.

    rp : 1D array
        Coordinate values for the control voluemes' "plus" interfaces.

    f : ND array
        Dependent variable values that correspond to the control volumes'
        center coordinates along ``axis``.

    axis : int, optional
        The axis index of ``f`` that corresponds to the ``x`` coordinate. The
        default is 0.

    Returns
    -------
    int_r : float | ND array
        The result of integration. If ``f`` is 1D, the returned value will be
        a scalar. Otherwise, an array with one fewer dimensions than ``f`` will
        be returned.
    """

    import numpy as np

    r = 0.5 * (rm + rp)

    new_axes = [*range(f.ndim)]
    new_axes.remove(axis)

    r = np.expand_dims(r, axis=new_axes)
    rm = np.expand_dims(rm, axis=new_axes)
    rp = np.expand_dims(rp, axis=new_axes)

    int_r = np.sum(4 * np.pi * r**2 * f * (rp - rm), axis=axis)

    return int_r


def param_combinations(params: list[str],
                       values: list[_ndarray]) -> list[dict]:
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
