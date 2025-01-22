from __future__ import annotations

import inspect
from warnings import warn
from typing import TYPE_CHECKING

import numpy as np
from scipy import linalg
from sksundae import ida

if TYPE_CHECKING:
    from typing import Callable, Any
    from numpy.typing import ArrayLike


class IDASolver(ida.IDA):
    pass


class IDAResult(ida.IDAResult):
    pass


def bandwidth(resfn: Callable, t0: float, y0: ArrayLike, yp0: ArrayLike,
              userdata: Any = None, return_pattern: bool = False,
              ) -> tuple[int | np.ndarray]:
    """
    Determine half bandwiths.

    Given a residual function `resfn`, this routine perturbs values in `y0`
    and `yp0` to approximate a Jacobian pattern. The pattern is then used to
    find the lower and upper half bandwidths.

    Parameters
    ----------
    resfn : Callable
        Residuals function like `f(t, y, yp, res[, userdata])`.
    t0 : float
        Initial time at which `y0` and `yp0` are given.
    y0 : ArrayLike, shape(m,)
        Initial state.
    yp0 : ArrayLike, shape(m,)
        Initial state time derivatives.
    userdata : Any | None, optional
        Optional argument to pass to `resfn`. The default is None.
    return_pattern : bool, optional
        If True, returns the Jacobian pattern along with the half bandwidths.
        The default is False.

    Returns
    -------
    lband : int
        Lower half bandwidth.
    uband : int
        Upper half bandwidth.
    j_pattern : np.ndarray, shape(m, m)
        Jacobian pattern. Only returned if `return_pattern=True`.

    Warnings
    --------
    UserWarning
        'resfn' signature has 5 inputs, but 'userdata=None'.

    Raises
    ------
    ValueError
        'resfn' signature must have either 4 or 5 inputs.

    """

    # Wrap resfn for cases w/ and w/o userdata
    signature = inspect.signature(resfn)

    if len(signature.parameters) == 4:
        wrapper = lambda t, y, yp, res: resfn(t, y, yp, res)
    elif len(signature.parameters) == 5:
        if userdata is None:
            warn("'resfn' signature has 5 inputs, but 'userdata=None'.")
        wrapper = lambda t, y, yp, res: resfn(t, y, yp, res, userdata)
    else:
        raise ValueError("'resfn' signature must have either 4 or 5 inputs.")

    # Perturbed variables
    y = np.asarray(y0, copy=True)
    yp = np.asarray(yp0, copy=True)

    # Initial residuals
    res = np.zeros_like(y)
    res_0 = np.zeros_like(y)

    wrapper(t0, y, yp, res_0)

    rng = np.random.default_rng(seed=42)
    rand = rng.random(2)

    # Jacobian pattern
    def j_pattern(j):
        y_store, yp_store = y[j], yp[j]

        y[j] += max(1e-6, 1e-6*y[j]) * rand[0]
        yp[j] += max(1e-6, 1e-6*yp[j]) * rand[1]
        wrapper(t0, y, yp, res)

        y[j], yp[j] = y_store, yp_store

        return np.where(res_0 - res != 0, 1, 0)

    j_cols = [j_pattern(j) for j in range(y.size)]
    j_pat = np.column_stack(j_cols)

    # Find lband and uband
    output = linalg.bandwidth(j_pat)

    if return_pattern:
        output += (j_pat,)

    return output
