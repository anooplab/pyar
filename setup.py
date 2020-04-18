#!/usr/bin/env python3
from distutils.core import setup

setup(name='PyAR',
      version='0.2',
      packages=['pyar', 'pyar.afir', 'pyar.data_analysis', 'pyar.interface', 'pyar.data'],
      scripts=['pyar/scripts/pyar-cli', 'pyar/scripts/pyar-optimiser', 'pyar/scripts/pyar-tabu',
               'pyar/scripts/pyar-clustering'],
      url='https://github.com/anooplab/pyar',
      license='GPl v3',
      author='Anoop et al',
      author_email='anoop@chem.iitkgp.ac.in',
      description='',
      requires=['numpy', 'sklearn', 'scipy', 'pandas', 'matplotlib'],
      )
