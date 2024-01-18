<picture>
  <source media="(prefers-color-scheme: dark)" srcset="./logos/dark.svg">
  <img alt="BatMods-lite logo" src="./logos/light.svg">
</picture> 

[![CI badge][cibadge]][cilink] &nbsp; [![style badge][pep8badge]][pep8link]

[cibadge]: https://github.com/NREL/BatMods-lite/actions/workflows/ci.yaml/badge.svg
[cilink]: https://github.com/NREL/BatMods-lite/actions/workflows/ci.yaml

[pep8badge]: https://img.shields.io/badge/code%20style-pep8-orange.svg
[pep8link]: https://www.python.org/dev/peps/pep-0008

BatMods-lite is a Python package that includes:

1) Library of pre-built battery models
2) Kinetic/transport properties for common battery materials

## Installing
We recommend using [Anaconda](https://www.anaconda.com/download) to create and manage your Python virtual environments. The following directions assume you are using Anaconda. 

After installing Anaconda and downloading the BatMods-lite repo files, open Anaconda Prompt (Windows) or Terminal (MacOS) and run the following to: (1) create an environment named "batmods," (2) activate your new environment, (3) install scikits-odes, and (4) install batmods-lite. 

1) ``conda create -n batmods python=3.10``
2) ``conda activate batmods``
3) ``conda install scikits.odes``
4) ``pip install ./batmods-lite`` or ``pip install -e ./batmods-lite``

**Notes:**
* Step (3) installs scikits.odes separately through ``conda install`` rather than ``pip install`` to simplify setup. Installing scikits.odes via ``pip`` complicates this step by requiring a pre-installed SUNDIALS and C++/Fortran compilers. 

* The ``-e`` flag in step (4) installs BatMods-lite in "editable" mode. You should use this if you plan to make changes to the package.

## Get Started
You are now ready to start running models. Run the following from your favorite terminal/IDE to see helpful documentation, examples, and more:

```python
import bmlite as bm

bm.docs()
```

## Formatting
BatMods-lite code follows the PEP8 style guide. It is the autors' preference to not adopt the more opinionated `black` formatting style. Please keep this in mind if you plan to contribute.

## Acknowledgements
This work was authored by the National Renewable Energy Laboratory (NREL), operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by funding from DOE's Vehicle Technologies Office (VTO) and Advanced Scientific Computing Research (ASCR) program. The research was performed using computational resources sponsored by the Department of Energy's Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory. The views expressed in the repository do not necessarily represent the views of the DOE or the U.S. Government.

## References
