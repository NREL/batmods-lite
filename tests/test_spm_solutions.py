import warnings

import pytest
import numpy as np
import bmlite as bm
import matplotlib.pyplot as plt


@pytest.fixture(scope='function')
def soln():

    warnings.filterwarnings('ignore')

    model = bm.SPM.Simulation()

    expr = bm.Experiment()
    expr.add_step('current_C', -2., (3600., 1.), limits=('voltage_V', 3.))
    expr.add_step('current_A', 0., (3600., 1.))

    soln = model.run(expr)

    return soln


def test_step_solution(soln):

    # solvetime works
    step_soln = soln.get_steps(0)
    assert step_soln.solvetime

    # bad plot
    with pytest.raises(KeyError):
        step_soln.simple_plot('fake', 'plot')

    # good plot
    with plt.ioff():
        step_soln.simple_plot('time_h', 'voltage_V')
        plt.close('all')

    # complex plots
    args = ('potentials', 'intercalation', 'pixels')
    with plt.ioff():
        step_soln.complex_plot(*args)
        plt.close('all')


def test_cycle_solution(soln):

    # solvetime works and times stacked correctly
    cycle_soln = soln.get_steps((0, 1))
    assert cycle_soln.solvetime
    assert all(np.diff(cycle_soln.t) >= 0.)

    # bad plot
    with pytest.raises(KeyError):
        cycle_soln.simple_plot('fake', 'plot')

    # good plot
    with plt.ioff():
        cycle_soln.simple_plot('time_h', 'voltage_V')
        plt.close('all')

    # complex plots
    args = ('potentials', 'intercalation', 'pixels')
    with plt.ioff():
        cycle_soln.complex_plot(*args)
        plt.close('all')
