import pytest

import bmlite as bm


@pytest.fixture(scope='session')
def sim():
    sim = bm.P2D.Simulation()
    return sim


@pytest.fixture(scope='session')
def sol(sim):
    exp = {
        'V_ext': 3.8,
        't_min': 0.0,
        't_max': 1350.0,
        'Nt': 150
    }

    sol = sim.run_CV(exp)
    return sol


def test_classname(sol):
    assert sol.classname == 'CVSolution'


def test_run_CV(sol):
    assert sol.success


def test_verify(sol):
    assert True
