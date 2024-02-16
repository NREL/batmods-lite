import os
from setuptools import setup

from bmlite import __version__

here = os.path.dirname(__file__)
with open(here + '/requirements.txt') as f:
    install_requires = f.readlines()

setup(
      name='batmods-lite',
      version=__version__,
      description='Pre-built physics-based battery models',
      author='Corey R. Randall',
      package_dir={'bmlite': 'bmlite'},
      classifiers=['Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Programming Language :: Python :: 3.10'
                   ],
      python_requires='>=3.10',
      install_requires=install_requires
     )
