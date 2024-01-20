import pytest

import bmlite as bm


@pytest.fixture(scope='session')
def sim():
    sim = bm.SPM.Simulation()
    return sim


@pytest.fixture(scope='session')
def sol(sim):
    exp = {
        'C_rate': -2.0,
        't_min': 0.0,
        't_max': 1350.0,
        'Nt': 150
    }

    sol = sim.run_CC(exp)
    return sol


def test_classname(sol):
    assert sol.classname == 'CCSolution'


def test_run_CC(sol):
    assert sol.success


def test_verify(sol):
    assert sol.verify(True)
