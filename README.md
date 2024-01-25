<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./images/dark.svg">
  <img alt="BatMods-lite logo" src="./images/light.svg">
</picture> 


<div>
[![CI badge][ci-badge]][ci-link] &nbsp; 
<img alt="tests badge" src="./images/tests.svg"> &nbsp; 
<img alt="coverage badge" src="./images/coverage.svg"> &nbsp; 
[![pep8 badge][pep-badge]][pep-link]
</div>

[ci-badge]: https://github.com/NREL/BatMods-lite/actions/workflows/ci.yaml/badge.svg
[ci-link]: https://github.com/NREL/BatMods-lite/actions/workflows/ci.yaml

[pep-badge]: https://img.shields.io/badge/code%20style-pep8-orange.svg
[pep-link]: https://www.python.org/dev/peps/pep-0008

</br>
BatMods-lite is a Python package that includes:

1) A library and API for pre-built battery models
2) Kinetic/transport properties for common battery materials

## Installing
We recommend using [Anaconda](https://www.anaconda.com/download) to create and manage your Python virtual environments. The following directions assume you are using Anaconda. 

After installing Anaconda and downloading or cloning the BatMods-lite repo files onto your local machine, open Anaconda Prompt (Windows) or Terminal (MacOS/Linux) and run the following to: (1) create an environment named "batmods," (2) activate your new environment, (3) install scikits.odes, and (4) install batmods-lite. 

1) ``conda create -n batmods python=3.10``
2) ``conda activate batmods``
3) ``conda install scikits.odes``
4) ``pip install ./batmods-lite`` or ``pip install -e ./batmods-lite``

**Notes:**
* Feel free to replace "batmods" in steps (1) and (2) with your preferred environment name.

* Step (3) uses ``conda install`` rather than ``pip install`` to simplify setup. Installing scikits.odes via ``pip`` requires a pre-installed SUNDIALS and C++/Fortran compilers. If you'd like to try to skip using ``conda``, please refer to the [scikits.odes documenation](https://scikits-odes.readthedocs.io/en/latest/).

* The ``-e`` flag in step (4) installs BatMods-lite in "editable" mode. You should use this if you plan to make changes to the package.

* The path ``./batmods-lite`` in step (4) assumes your terminal's working directory is the parent directory for the "batmods-lite" folder. If this is not true, please ``cd`` into the correct directory, or change the path to point to the correct folder.

## Get Started
You are now ready to start running models. Run the following from your favorite terminal/IDE to see helpful documentation, examples, and more:

```python
import bmlite as bm

bm.docs()
```

**Notes:**
* If you are new to Python, check out [Spyder IDE](https://www.spyder-ide.org/). Spyder is a powerful interactive development environment. BatMods-lite is programmed almost entirely using the Spyder IDE. 

* If you choose to use Spyder, you can install it using ``conda install spyder``. Afterward, you can find the application on your system, or simply type ``spyder`` in Anaconda Prompt or Terminal to open an instance.

## Formatting
With minor exceptions, BatMods-lite code follows the [PEP8 style guide](https://www.python.org/dev/peps/pep-0008). Specifically, adding extra spaces around parentheses, and under- or over-indenting multi-line expressions is allowed when it improves readability or avoids the 80-character line limit. 

Please be aware that it is the authors' preference to not adopt the more opinionated [black formatting style](https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html). Make sure that you do not autoformat any files in this style if you choose to contribute.

## Acknowledgements
This work was authored by the National Renewable Energy Laboratory (NREL), operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by funding from DOE's Vehicle Technologies Office (VTO) and Advanced Scientific Computing Research (ASCR) program. The research was performed using computational resources sponsored by the Department of Energy's Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory. The views expressed in the repository do not necessarily represent the views of the DOE or the U.S. Government.
