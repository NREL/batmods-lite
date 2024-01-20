<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./images/dark.svg">
  <img alt="BatMods-lite logo" src="./images/light.svg">
</picture> 

<div>
<a href="https://github.com/NREL/BatMods-lite/actions/workflows/ci.yaml">
  <img alt="CI badge" src="https://github.com/NREL/BatMods-lite/actions/workflows/ci.yaml/badge.svg">
</a> &nbsp; 
<img alt="coverage badge" src="./images/coverage.svg"> &nbsp; 
<a href="https://www.python.org/dev/peps/pep-0008">
  <img alt="style badge" src="https://img.shields.io/badge/code%20style-pep8-orange.svg">
</a>
</div>

BatMods-lite is a Python package that includes:

1) A library of pre-built battery models
2) Kinetic/transport properties for common battery materials

## Installing
We recommend using [Anaconda](https://www.anaconda.com/download) to create and manage your Python virtual environments. The following directions assume you are using Anaconda. 

After installing Anaconda and downloading the BatMods-lite repo files, open Anaconda Prompt (Windows) or Terminal (MacOS/Linux) and run the following to: (1) create an environment named "batmods," (2) activate your new environment, (3) install scikits-odes, and (4) install batmods-lite. 

1) ``conda create -n batmods python=3.10``
2) ``conda activate batmods``
3) ``conda install scikits.odes``
4) ``pip install ./batmods-lite`` or ``pip install -e ./batmods-lite``

**Notes:**
* Feel free to replace "batmods" in steps (1) and (2) with your preferred environment name.

* Step (3) uses ``conda install`` rather than ``pip install`` to simplify the setup. Installing scikits.odes via ``pip`` requires a pre-installed SUNDIALS and C++/Fortran compilers. 

* The ``-e`` flag in step (4) installs BatMods-lite in "editable" mode. You should use this if you plan to make changes to the package.

## Get Started
You are now ready to start running models. Run the following from your favorite terminal/IDE to see helpful documentation, examples, and more:

```python
import bmlite as bm

bm.docs()
```

**Notes:**
* If you are new to Python, check out [Spyder IDE](https://www.spyder-ide.org/). Spyder is a powerful interactive development environment. BatMods-lite is built almost entirely using the Spyder IDE. 

* If you choose to use Spyder, you can install it using ``conda install spyder``. Afterward, you can find the application on your system, or simply type ``spyder`` in Anaconda Prompt or Terminal to open an instance.

## Formatting
With minor exceptions, BatMods-lite code follows the [PEP8 style guide](https://www.python.org/dev/peps/pep-0008). Specifically, adding extra spaces around brackets/parentheses, and under- or over-indenting multi-line expressions is allowed when it improves readability. 

Please be aware that it is the authors' preference to not adopt the more opinionated [black formatting style](https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html). Make sure that you do not autoformat any files in this style if you choose to contribute.

## Acknowledgements
This work was authored by the National Renewable Energy Laboratory (NREL), operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by funding from DOE's Vehicle Technologies Office (VTO) and Advanced Scientific Computing Research (ASCR) program. The research was performed using computational resources sponsored by the Department of Energy's Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory. The views expressed in the repository do not necessarily represent the views of the DOE or the U.S. Government.
