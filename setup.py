from setuptools import setup

setup(
      name='batmods-lite',
      version='0.0.1',
      description='Pre-built physics-based battery models',
      author='Corey R. Randall',
      package_dir={'bmlite': 'bmlite'},
      classifiers=['Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Programming Language :: Python :: 3.10'
                   ],
      python_requires='>=3.10',
      install_requires=['numpy', 'pandas', 'openpyxl', 'scipy', 'matplotlib',
                        'ruamel.yaml']
     )
