import numpy as np
import bmlite as bm


def test_grad_x():

    x = np.linspace(0, 10, 11)
    f = np.sin(x)

    grad_x = bm.mathutils.grad_x(x, f)

    df_dx = (f[1:] - f[:-1]) / (x[1:] - x[:-1])

    assert np.allclose(grad_x, df_dx)


def test_grad_r():

    r = np.linspace(0, 10, 11)
    f = np.sin(r)

    grad_r = bm.mathutils.grad_r(r, f)

    df_dr = (f[1:] - f[:-1]) / (r[1:] - r[:-1])

    assert np.allclose(grad_r, df_dr)


def test_div_x():

    xm = np.linspace(0, 10, 11)
    xp = np.linspace(1, 11, 11)

    f = np.sin(np.hstack([xm, xp[-1]]))

    div_x = bm.mathutils.div_x(xm, xp, f)

    df_dx = (f[1:] - f[:-1]) / (xp - xm)

    assert np.allclose(div_x, df_dx)


def test_div_r():

    rm = np.linspace(0, 10, 11)
    rp = np.linspace(1, 11, 11)
    r = 0.5 * (rm + rp)

    f = np.sin(np.hstack([rm, rp[-1]]))

    div_r = bm.mathutils.div_r(rm, rp, f)

    df_dr = 1 / r**2 * (rp**2 * f[1:] - rm**2 * f[:-1]) / (rp - rm)

    assert np.allclose(div_r, df_dr)


def test_int_x():

    xm = np.linspace(0, 10, 11)
    xp = np.linspace(1, 11, 11)
    x = 0.5 * (xm + xp)

    f = np.sin(x)

    int_x = bm.mathutils.int_x(xm, xp, f)

    integral = np.sum(f * (xp - xm))

    assert int_x == integral


def test_int_r():

    rm = np.linspace(0, 10, 11)
    rp = np.linspace(1, 11, 11)
    r = 0.5 * (rm + rp)

    f = np.sin(r)

    int_r = bm.mathutils.int_r(rm, rp, f)

    integral = np.sum(4 * np.pi * r**2 * f * (rp - rm))

    assert int_r == integral


def test_param_combinations():

    params = ['a', 'b']
    values = [np.linspace(0, 1, 2), np.linspace(3, 4, 2)]

    params_list = bm.mathutils.param_combinations(params, values)

    check_list = []
    for a in values[0]:
        new_combination = {}
        new_combination['a'] = a
        for b in values[1]:
            new_combination['b'] = b
            check_list.append(new_combination.copy())

    assert params_list == check_list
