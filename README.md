<picture>
  <source media='(prefers-color-scheme: dark)' srcset='./images/dark.svg'
  style='width: 75%; min-width: 250px; max-width: 500px;'>
  <img alt='BATMODS-lite logo' src='./images/light.svg'
  style='width: 75%; min-width: 250px; max-width: 500px;'>
</picture> 

[![CI][ci-b]][ci-l] &nbsp;
![tests][test-b] &nbsp;
![coverage][cov-b] &nbsp;
[![pep8][pep-b]][pep-l]

[ci-b]: https://github.com/NREL/BATMODS-lite/actions/workflows/ci.yaml/badge.svg
[ci-l]: https://github.com/NREL/BATMODS-lite/actions/workflows/ci.yaml

[test-b]: ./images/tests.svg
[cov-b]: ./images/coverage.svg

[pep-b]: https://img.shields.io/badge/code%20style-pep8-orange.svg
[pep-l]: https://www.python.org/dev/peps/pep-0008

## Summary
Battery Analysis and Training Models for Optimization and Degradation Studies (BATMODS) is a Python package with an API for pre-built battery models. The original purpose of the package was to quickly generate synthetic data for machine learning models to train with. However, the models are generally useful for any battery simulations or analysis. BATMODS-lite includes the following: 

1) A library and API for pre-built battery models
2) Kinetic/transport properties for common battery materials

## Installation
BATMODS-lite is only available via GitHub. Please clone the repo or download the files. Also, make sure you have a Python installation with a version >=3.10. If you are new to Python, we recommend using [Anaconda](https://www.anaconda.com/download) to set up your installation. 

Once the files are available on your machine, use your terminal to navigate into the folder and execute one of the following depending on your installation preference.

```
pip install .             (basic installation)
pip install -e .[dev]     (editable installation with developer options)
```

The editable installation is useful if you plan to make changes to your local package. It ensures that any changes are immediately available each time the package is imported, without needing to reinstall. The developer options will likely be helpful if you are modifying the package. For more info on pre-build developer routines look at `noxfile.py`.

## Get Started
The API is organized around three main classes that allow you to construct simulations, define experiments, and interact with solutions. Two basic examples are given below. These demonstrate a 2C discharge for both the single particle model (SPM) and pseudo-2D (P2D) model. Note that the experiment class interfaces with all simulations. The simulations and their respective solutions, however, will depend on the model subpackage they are loaded from.

```python
# Single particle model example
import bmlite as bm

sim = bm.SPM.Simulation()

expr = bm.Experiment()
expr.add_step('current_C', 2., (1350., 150))

soln = sim.run(expr)
soln.simple_plot('time_s', 'voltage_V')
```

```python
# Pseudo-2D model example
import bmlite as bm

sim = bm.P2D.Simulation()

expr = bm.Experiment()
expr.add_step('current_C', 2., (1350., 150))

soln = sim.run(expr)
soln.simple_plot('time_s', 'voltage_V')
```

**Notes:**
* If you are new to Python, check out [Spyder IDE](https://www.spyder-ide.org/). Spyder is a powerful interactive development environment (IDE) that can make programming in Python more approachable to new users.

## Citing this Work

## Acknowledgements

## Contributing
If you plan to contribute, please use `nox` to ensure the package is developed and tested in a consistent way between local copies and continuous integration (CI) runs. You can read through sessions in `noxfile.py` to see what all functionality is already set up for you. In a brief overview, the following shows how to lint, spellcheck, and run tests.

```
nox -s linter
nox -s codespell
nox -s tests
```

All of these and more are run with `nox -s pre-commit` which should be used before all commits back to the remote repository. Running the pre-commit session locally saves time by helping to catch errors that may occur during the CI builds and tests.

## Disclaimer
This work was authored by the National Renewable Energy Laboratory (NREL), operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE). The views expressed in the repository do not necessarily represent the views of the DOE or the U.S. Government.
