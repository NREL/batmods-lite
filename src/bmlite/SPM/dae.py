"""
DAE Module
----------
This module includes the system of differential algebraic equations (DAE) for
the SPM model. In addition, the ``bandwidth`` function is included in this
module, which helps determine the lower and upper bandwidths of ``residuals``
so the ``'band'`` linear solver option can be used in the ``IDASolver`` class.

"""

import numpy as np

if not hasattr(np, 'concat'):  # pragma: no cover
    np.concat = np.concatenate


def sign(x):
    """Return +1 for x >= 0, -1 for x < 0."""
    return (x >= 0).astype(float) * 2 - 1


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
    outputs : tuple[np.ndarray]
        If the experimental step `mode` is set to `post`, then the following
        post-processed variables will be returned in a tuple. Otherwise,
        returns None.

        ========= =================================================
        Variable  Description [units] (*type*)
        ========= =================================================
        sdot_an   anode Li+ production [kmol/m^3/s] (*1D array*)
        sdot_ca   cathode Li+ production [kmol/m^3/s] (*1D array*)
        ========= =================================================

    """

    from ..mathutils import grad_r, div_r

    # Break inputs into separate objects
    sim, exp = inputs

    c, bat, el, an, ca = sim.c, sim.bat, sim.el, sim.an, sim.ca

    # Simulation temperature
    T = bat.temp

    # Organize values from sv
    phi_an = sv[an.ptr['phi_ed']]
    phi_el = sv[el.ptr['phi_el']]
    phi_ca = sv[ca.ptr['phi_ed']]

    if 'Hysteresis' in an._submodels:
        hyst_an = sv[an.ptr['hyst']]
        Hyst_an = an.M_hyst*hyst_an
    else:
        Hyst_an = 0.

    if 'Hysteresis' in ca._submodels:
        hyst_ca = sv[ca.ptr['hyst']]
        Hyst_ca = ca.M_hyst*hyst_ca
    else:
        Hyst_ca = 0.

    xs_an = sv[an.r_ptr['Li_ed']]
    xs_ca = np.flip(sv[ca.r_ptr['Li_ed']])

    Li_an = xs_an*an.Li_max
    Li_ca = xs_ca*ca.Li_max

    # Anode -------------------------------------------------------------------

    # Reaction current
    eta = phi_an - phi_el - (an.get_Eeq(xs_an[-1]) + Hyst_an)
    fluxdir_an = -np.sign(eta)

    i0 = an.get_i0(xs_an[-1], el.Li_0, T, fluxdir_an)
    sdot_an = i0 / c.F * (  np.exp( an.alpha_a*c.F*eta / c.R / T)
                          - np.exp(-an.alpha_c*c.F*eta / c.R / T)  )

    # Weighted solid particle properties
    Ds_an = an._wtm*an.get_Ds(xs_an[:-1], T, fluxdir_an) \
          + an._wtp*an.get_Ds(xs_an[1:], T, fluxdir_an)

    # Solid-phase COM (differential)
    Js_an = np.concat([[0.], Ds_an*grad_r(an.r, Li_an), [-sdot_an]])

    res[an.r_ptr['Li_ed']] = an.Li_max*svdot[an.r_ptr['Li_ed']] \
                           - div_r(an.rm, an.rp, Js_an)

    # Solid-phase COC (algebraic)
    res[an.ptr['phi_ed']] = phi_an - 0.

    # Hysteresis (differential)
    if 'Hysteresis' in an._submodels:
        res[an.ptr['hyst']] = svdot[an.ptr['hyst']] \
            - np.abs(sdot_an*c.F*an.g_hyst / 3600. / bat.cap) \
            * (sign(sdot_an) - hyst_an)

    # Cathode -----------------------------------------------------------------

    # Reaction current
    eta = phi_ca - phi_el - (ca.get_Eeq(xs_ca[-1]) + Hyst_ca)
    fluxdir_ca = -np.sign(eta)

    i0 = ca.get_i0(xs_ca[-1], el.Li_0, T, fluxdir_ca)
    sdot_ca = i0 / c.F * (  np.exp( ca.alpha_a*c.F*eta / c.R / T)
                          - np.exp(-ca.alpha_c*c.F*eta / c.R / T)  )

    # Weighted solid particle properties
    Ds_ca = ca._wtm*ca.get_Ds(xs_ca[:-1], T, fluxdir_ca) \
          + ca._wtp*ca.get_Ds(xs_ca[1:], T, fluxdir_ca)

    # Solid-phase COM (differential)
    Js_ca = np.concat([[0.], Ds_ca*grad_r(ca.r, Li_ca), [-sdot_ca]])

    res[ca.r_ptr['Li_ed']] = ca.Li_max*svdot[ca.r_ptr['Li_ed']] \
                           - np.flip(div_r(ca.rm, ca.rp, Js_ca))

    # Hysteresis (differential)
    if 'Hysteresis' in ca._submodels:
        res[ca.ptr['hyst']] = svdot[ca.ptr['hyst']] \
            - np.abs(sdot_ca*c.F*ca.g_hyst / 3600. / bat.cap) \
            * (sign(sdot_ca) - hyst_ca)

    # External current [A/m^2]
    i_ext = sdot_an*an.A_s*an.thick*c.F

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
                              + value(t) / bat.area
        res[el.ptr['phi_el']] = sdot_an*an.A_s*an.thick*c.F \
                              - value(t) / bat.area

    elif mode == 'current' and units == 'C':
        res[ca.ptr['phi_ed']] = sdot_ca*ca.A_s*ca.thick*c.F \
                              + value(t)*bat.cap / bat.area
        res[el.ptr['phi_el']] = sdot_an*an.A_s*an.thick*c.F \
                              - value(t)*bat.cap / bat.area

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
    if mode == 'post':
        return sdot_an, sdot_ca
