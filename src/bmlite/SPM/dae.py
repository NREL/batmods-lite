"""
DAE Module
----------
This module includes the system of differential algebraic equations (DAE) for
the SPM model. In addition, the ``bandwidth`` function is included in this
module, which helps determine the lower and upper bandwidths of ``residuals``
so the ``'band'`` linear solver option can be used in the ``IDASolver`` class.

"""

import numpy as np


def bandwidth(sim: object) -> tuple[int | np.ndarray]:
    """
    Determine the DAE system's bandwidth and Jacobian pattern.

    Numerically determines the bandwidth and Jacobian pattern of the residual
    function by perturbating each ``y`` and ``yp`` term and determining which
    ``dres/dy`` and ``dres/dy'`` terms are non-zero. The bandwidth is required
    to use the "band" option for IDA's linear solver, which speeds up each
    integration step compared to the "dense" linear solver.

    Parameters
    ----------
    inputs : SPM Simulation object
        An instance of the SPM model simulation. See
        :class:`bmlite.SPM.Simulation`.

    Returns
    -------
    lband : int
        Lower bandwidth from the residual function's Jacobian pattern.
    uband : int
        Upper bandwidth from the residual function's Jacobian pattern.
    j_pat : 2D array
        Residual function Jacobian pattern, as an array of ones and zeros.

    """

    # Jacobian size
    N = sim._sv0.size

    # Fake OCV experiment
    expr = {
        'mode': 'current',
        'units': 'C',
        'value': lambda t: 0.,
    }

    # Perturbed variables
    jac = np.zeros([N, N])
    sv = sim._sv0.copy()
    svdot = sim._svdot0.copy()
    res = np.zeros_like(sv)

    residuals(0, sv, svdot, res, (sim, expr))
    res_0 = res.copy()

    for j in range(N):
        dsv = np.copy(sv)
        res = np.copy(res)
        svdot = np.copy(svdot)

        dsv[j] = sv[j] + max(1e-6, 1e-6*sv[j])
        residuals(0, dsv, svdot, res, (sim, expr))
        dres = res_0 - res

        jac[:, j] = dres

    for j in range(N):
        sv = np.copy(sv)
        res = np.copy(res)
        dsvdot = np.copy(svdot)

        dsvdot[j] = svdot[j] + max(1e-6, 1e-6*svdot[j])
        residuals(0., sv, dsvdot, res, (sim, expr))
        dres = res_0 - res

        jac[:, j] += dres

    # Find lband and uband
    lband = 0
    uband = 0

    for i in range(jac.shape[0]):

        l_inds = np.where(abs(jac[i, :i]) > 0)[0]
        if len(l_inds) >= 1 and i - l_inds[0] > lband:
            lband = int(i - l_inds[0])

        u_inds = i + np.where(abs(jac[i, i:]) > 0)[0]
        if len(u_inds) >= 1 and u_inds[-1] - i > uband:
            uband = int(u_inds[-1] - i)

    # Make Jacobian pattern of zeros and ones
    j_pat = np.zeros_like(jac)
    j_pat[jac != 0] = 1

    return lband, uband, j_pat


def residuals(t: float, sv: np.ndarray, svdot: np.ndarray, res: np.ndarray,
              inputs: tuple[object, dict]) -> None | tuple[np.ndarray]:
    """
    The DAE residuals ``res = M*y' - f(t, y)`` for the SPM model.

    Parameters
    ----------
    t : float
        Value of time [s].
    sv : 1D array
        Solution/state variables at time ``t``.
    svdot : 1D array
        Solution/state variable time derivatives at time ``t``.
    res : 1D array
        An array the same size as ``sv`` and ``svdot``. The values are filled
        in with ``res = M*y' - f(t, y)`` inside this function.
    inputs : (sim : SPM Simulation object, exp : experiment dict)
        The simulation object and experimental details dictionary inputs that
        describe the specific battery and experiment to simulate.

    Returns
    -------
    None
        If no ``sim._flags`` are ``True``.
    res : 1D array
        Array of residuals if ``sim._flags['band'] = True``.
    outputs : tuple[1D array]
        If ``sim._flags['post'] = True`` then ``outputs`` is returned, which
        includes post-processed values. These can help verify the governing
        equations and boundary conditions are satisfied. They can also be
        useful for interpreting causes of good/bad battery performance. The
        order and description of the arrays is given below:

        ========== =======================================================
        Variable   Description [units] (*type*)
        ========== =======================================================
        res        residuals ``res = M*y' - f(t, y)`` [units] (*1D array*)
        sdot_an    anode Li+ production rate [kmol/m^3/s] (*float*)
        sdot_ca    cathode Li+ production rate [kmol/m^3/s] (*float*)
        ========== =======================================================

    """

    from .. import Constants
    from ..math import grad_r, div_r

    c = Constants()

    # Break inputs into separate objects
    sim, exp = inputs

    bat, el, an, ca = sim.bat, sim.el, sim.an, sim.ca

    # Simulation temperature
    T = bat.temp

    # Organize values from sv
    phi_an = sv[an.ptr['phi_ed']]
    phi_el = sv[el.ptr['phi_el']]
    phi_ca = sv[ca.ptr['phi_ed']]

    xs_an = sv[an.r_ptr['Li_ed']]
    xs_ca = sv[ca.r_ptr['Li_ed']]

    Li_an = xs_an*an.Li_max
    Li_ca = xs_ca*ca.Li_max

    # Anode -------------------------------------------------------------------

    # Reaction current
    eta = phi_an - phi_el - an.get_Eeq(xs_an[-1], T)

    i0 = an.get_i0(xs_an[-1], el.Li_0, T)
    sdot_an = i0 / c.F * (  np.exp( an.alpha_a * c.F * eta / c.R / T)
                          - np.exp(-an.alpha_c * c.F * eta / c.R / T)  )

    # Weighted solid particle properties
    wt_m = 0.5*(an.rp[:-1] - an.rm[:-1]) / (an.r[1:] - an.r[:-1])
    wt_p = 0.5*(an.rp[1:] - an.rm[1:]) / (an.r[1:] - an.r[:-1])

    Ds_an = wt_m*an.get_Ds(xs_an[:-1], T) + wt_p*an.get_Ds(xs_an[1:], T)

    # Solid-phase COM (differential)
    Js_an = np.concat([[0.], Ds_an*grad_r(an.r, Li_an), [-sdot_an]])

    res[an.r_ptr['Li_ed']] = an.Li_max*svdot[an.r_ptr['Li_ed']] \
                           - div_r(an.rm, an.rp, Js_an)

    # Solid-phase COC (algebraic)
    res[an.ptr['phi_ed']] = phi_an - 0.

    # Cathode -----------------------------------------------------------------

    # Reaction current
    eta = phi_ca - phi_el - ca.get_Eeq(xs_ca[-1], T)

    i0 = ca.get_i0(xs_ca[-1], el.Li_0, T)
    sdot_ca = i0 / c.F * (  np.exp( ca.alpha_a*c.F*eta / c.R / T)
                          - np.exp(-ca.alpha_c*c.F*eta / c.R / T)  )

    # Weighted solid particle properties
    wt_m = 0.5*(ca.rp[:-1] - ca.rm[:-1]) / (ca.r[1:] - ca.r[:-1])
    wt_p = 0.5*(ca.rp[1:] - ca.rm[1:]) / (ca.r[1:] - ca.r[:-1])

    Ds_ca = wt_m*ca.get_Ds(xs_ca[:-1], T) + wt_p*ca.get_Ds(xs_ca[1:], T)

    # Solid-phase COM (differential)
    Js_ca = np.concat([[0.], Ds_ca*grad_r(ca.r, Li_ca), [-sdot_ca]])

    res[ca.r_ptr['Li_ed']] = ca.Li_max*svdot[ca.r_ptr['Li_ed']] \
                           - div_r(ca.rm, ca.rp, Js_ca)

    # External current [A/m^2]
    i_ext = -sdot_an*an.A_s*an.thick*c.F

    # Boundary conditions -----------------------------------------------------
    mode = exp['mode']
    units = exp['units']
    value = exp['value']

    voltage_V = phi_ca
    current_A = i_ext*bat.area
    power_W = current_A*voltage_V

    # Cathode - Solid-phase COC (algebraic)
    # Electrolyte - potential (algebraic)
    if mode == 'current' and units == 'A':
        res[ca.ptr['phi_ed']] = sdot_ca*ca.A_s*ca.thick*c.F \
                              - value(t) / bat.area
        res[el.ptr['phi_el']] = sdot_an*an.A_s*an.thick*c.F \
                              + value(t) / bat.area

    elif mode == 'current' and units == 'C':
        res[ca.ptr['phi_ed']] = sdot_ca*ca.A_s*ca.thick*c.F \
                              - value(t)*bat.cap / bat.area
        res[el.ptr['phi_el']] = sdot_an*an.A_s*an.thick*c.F \
                              + value(t)*bat.cap / bat.area

    elif mode == 'voltage':
        res[ca.ptr['phi_ed']] = voltage_V - value(t)
        res[el.ptr['phi_el']] = sdot_an*an.A_s*an.thick \
                              + sdot_ca*ca.A_s*ca.thick

    elif mode == 'power':
        res[ca.ptr['phi_ed']] = power_W - value(t)
        res[el.ptr['phi_el']] = sdot_an*an.A_s*an.thick \
                              + sdot_ca*ca.A_s*ca.thick

    # Events tracking ---------------------------------------------------------
    total_time = sim._t0 + t

    exp['events'] = {
        'time_s': total_time,
        'time_min': total_time / 60.,
        'time_h': total_time / 3600.,
        'current_A': current_A,
        'current_C': current_A / bat.cap,
        'voltage_V': voltage_V,
        'power_W': power_W,
    }

    # Returns -----------------------------------------------------------------
    if sim._flags['post']:
        return sdot_an, sdot_ca
