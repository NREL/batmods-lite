<!-- <img alt='Logo' style='width: 75%; min-width: 250px; max-width: 500px;' 
 src='https://github.com/NREL/batmods-lite/blob/main/images/dark.png?raw=true#gh-dark-mode-only'/>
<img alt='Logo' style='width: 75%; min-width: 250px; max-width: 500px;' 
 src='https://github.com/NREL/batmods-lite/blob/main/images/light.png?raw=true#gh-light-mode-only'/> -->

# BATMODS-lite

[![CI][ci-b]][ci-l] &nbsp;
![tests][test-b] &nbsp;
![coverage][cov-b] &nbsp;
[![pep8][pep-b]][pep-l]

[ci-b]: https://github.com/NREL/batmods-lite/actions/workflows/ci.yml/badge.svg
[ci-l]: https://github.com/NREL/batmods-lite/actions/workflows/ci.yml

[test-b]: https://github.com/NREL/batmods-lite/blob/main/images/tests.svg?raw=true
[cov-b]: https://github.com/NREL/batmods-lite/blob/main/images/coverage.svg?raw=true

[pep-b]: https://img.shields.io/badge/code%20style-pep8-orange.svg
[pep-l]: https://www.python.org/dev/peps/pep-0008

## Summary
Battery Analysis and Training Models for Optimization and Design Studies (BATMODS) is a Python package with an API for pre-built battery models. The original purpose of the package was to quickly generate synthetic data for machine learning models to train with. However, the models are generally useful for any battery simulations or analysis. BATMODS-lite includes the following: 

1) A library and API for pre-built battery models
2) Kinetic/transport properties for common battery materials

## Installation
BATMODS-lite is only available via GitHub. Please clone the repo or download the files. Also, make sure you have a Python installation with a version >=3.10. If you are new to Python, we recommend using [Anaconda](https://www.anaconda.com/download) to set up your installation. 

Once the files are available on your machine, use your terminal to navigate into the folder and execute one of the following depending on your installation preference.

```
pip install .             (basic installation)
pip install -e .[dev]     (editable installation with developer options)
```

The editable installation is useful if you plan to make changes to your local package. It ensures that any changes are immediately available each time the package is imported, without needing to reinstall. The developer options will likely be helpful if you are modifying the package.

## Get Started
The API is organized around three main classes that allow you to construct simulations, define experiments, and interact with solutions. Two basic examples are given below. These demonstrate a 2C discharge for both the single particle model (SPM) and pseudo-2D (P2D) model. Note that the experiment class interfaces with all simulations. The simulations and their respective solutions, however, will depend on the model subpackage they are loaded from. For a more detailed tutorial, please check `docs/source/examples`. If you installed the editable version with developer options, you can also build the documentation locally using `nox -s docs`.

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
* Another friendly option for getting started in Python is to use [Jupyter Notebooks](https://jupyter.org/). We write our examples in Jupyter Notebooks since they support both markdown blocks for explanations and executable code blocks.
* Python, Spyder, and Jupyter Notebooks can be setup using [Anaconda](https://www.anaconda.com/download/success). Anaconda provides a convenient way for new users to get started with Python due to its friendly graphical installer and environment manager.

## Citing this Work
This work was authored by researchers at the National Laboratory of the Rockies (NLR). If you use this package in your work, please include the following citation:

> Randall, Corey R. "BATMODS-lite: Packaged battery models and material properties [SWR-25-108]." Computer software, Jun. 2025. url: [github.com/NREL/batmods-lite](https://github.com/NREL/batmods-lite). doi: [waiting for doi...](??).

For convenience, we also provide the following for your BibTex:

```
@misc{randall2025bmlite,
  author = {Randall, Corey R.},
  title = {{BATMODS-lite: Packaged battery models and material properties [SWR-25-108]}},
  url = {github.com/NREL/batmods-lite},
  month = {Jun.},
  year = {2025},
  doi = {waiting for doi...},
}
```

## Contributing
If you'd like to contribute to this package, please look through the existing [issues](https://github.com/NREL/batmods-lite/issues). If the bug you've caught or the feature you'd like to add isn't already being worked on, please submit a new issue before getting started. 

<!-- You should also read through the [developer guidelines](https://batmods-lite.readthedocs.io/latest/development). Need to get on readthedocs first... -->

## Disclaimer
This work was authored by the National Laboratory of the Rockies (NLR), operated by Alliance for Energy Innovation, LLC, for the U.S. Department of Energy (DOE). The views expressed in the repository do not necessarily represent the views of the DOE or the U.S. Government.
