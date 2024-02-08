from setuptools import setup
import os

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, "requirements.txt")) as f:
    install_requires = f.readlines()

with open(os.path.join(here, "bmlite", "version.py"), encoding="utf-8") as f:
    version = f.read()

version = version.split('=')[-1].strip().strip('"').strip("'")

setup(
      name='batmods-lite',
      version=version,
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
