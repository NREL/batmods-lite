import os, pytest
from pathlib import Path

import numpy as np
import bmlite as bm
from ruamel.yaml import YAML
import matplotlib.pyplot as plt


class Values(object):

    def __init__(self) -> None:
        self.t = None
        self.y = None
        self.ydot = None


class Roots(object):

    def __init__(self) -> None:
        self.t = None
        self.y = None
        self.ydot = None


class IDASol(object):

    def __init__(self) -> None:
        self.values = Values()
        self.values.t = np.linspace(0., 1350., 150)
        self.values.y = np.ones([10, 10])
        self.values.ydot = np.zeros([10, 10])
        self.roots = Roots()
        self.flag = 0
        self.onroot = False
        self.message = 'Successful solve.'


@pytest.fixture(scope='session')
def sim():
    sim = bm.P2D.Simulation()
    return sim


@pytest.fixture(scope='session')
def exp():
    directory = os.path.dirname(__file__) + '/../bmlite/P2D/default_exps/'
    path = Path(directory + 'constant_current.yaml')
    yaml = YAML()

    exp = yaml.load(path)
    return exp


@pytest.fixture(scope='session')
def sol(sim, exp):
    sol = bm.P2D.solutions.BaseSolution(sim, exp)
    return sol


@pytest.fixture(scope='session')
def idasol():
    idasol = IDASol()
    return idasol


@pytest.fixture(scope='session')
def dictsol():
    dictsol = {'t': np.linspace(0., 1350., 150),
               'y': np.ones([10, 10]),
               'ydot': np.zeros([10, 10]),
               'success': True,
               'onroot': False,
               'message': 'Successful solve.',
               'solvetime': 60.0
               }

    return dictsol


@pytest.fixture(scope='session')
def plotsol(sim):

    N = sim.sv_0.size
    M = sim.svdot_0.size

    dictsol = {'t': np.linspace(0., 1350., 150),
               'y': np.tile(sim.sv_0, 150).reshape(150, N),
               'ydot': np.tile(sim.svdot_0, 150).reshape(150, M),
               'success': True,
               'onroot': False,
               'message': 'Successful solve.',
               'solvetime': 60.0
               }

    return dictsol


def test_classname(sol):
    assert sol.classname == 'BaseSolution'


def test_ida_fill(sol, idasol):
    sol.ida_fill(idasol, 60.0)

    assert all(sol.t == idasol.values.t)
    assert np.all(sol.y == idasol.values.y)
    assert np.all(sol.ydot == idasol.values.ydot)
    assert sol.success == bool(idasol.flag >= 0)
    assert sol.onroot == bool(not isinstance(idasol.roots.t, type(None)))
    assert sol.message == idasol.message
    assert sol._solvetime == 60.0


def test_dict_fill(sol, dictsol):
    sol.dict_fill(dictsol)

    assert all(sol.t == dictsol['t'])
    assert np.all(sol.y == dictsol['y'])
    assert np.all(sol.ydot == dictsol['ydot'])
    assert sol.success == dictsol['success']
    assert sol.onroot == dictsol['onroot']
    assert sol.message == dictsol['message']
    assert sol._solvetime == dictsol['solvetime']


def test_solvetime(sol):
    sol._solvetime = 60.0

    assert sol.solvetime('s') == f'Solve time: {sol._solvetime:.3f} s'
    assert sol.solvetime('min') == f'Solve time: {sol._solvetime/60:.3f} min'
    assert sol.solvetime('h') == f'Solve time: {sol._solvetime/3600:.3f} h'


def test_report(sol):
    sol.report()
    assert True


def test_to_dict(sol, dictsol):
    sol.dict_fill(dictsol)
    testdict = sol.to_dict()

    assert testdict == dictsol


def test_plot(sol, plotsol):
    sol.dict_fill(plotsol)

    args = ('current', 'voltage', 'power', 'ivp', 'potentials', 'electrolyte',
            'intercalation', 'pixels')

    for a in args:
        assert callable(getattr(bm.P2D.postutils, a))
        sol.plot(a)

        plt.close('all')


# def test_slice_and_save(sol, plotsol):
#     sol.dict_fill(plotsol)

#     directory = os.path.dirname(__file__)
#     sol.slice_and_save(directory + '/scrap_data/test_p2d')
#     assert os.path.exists(directory + '/scrap_data/test_p2d.npz')

#     data = np.load(directory + '/scrap_data/test_p2d.npz')
#     data.close()

#     with pytest.raises(Exception):
#         sol.slice_and_save(directory + '/scrap_data/test_p2d')

#     sol.slice_and_save(directory + '/scrap_data/test_p2d', overwrite=True)
#     os.remove(directory + '/scrap_data/test_p2d.npz')
