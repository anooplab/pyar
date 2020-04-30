#!/usr/bin/env python3
from distutils.core import setup

setup(name='PyAR',
      version='0.3',
      packages=['pyar', 'pyar.afir', 'pyar.data_analysis', 'pyar.interface', 'pyar.data'],
      scripts=['pyar/scripts/pyar-cli', 'pyar/scripts/pyar-optimiser', 'pyar/scripts/pyar-tabu',
               'pyar/scripts/pyar-clustering'],
      url='https://github.com/anooplab/pyar',
      license='GPl v3',
      author='Anoop et al',
      author_email='anoop@chem.iitkgp.ac.in',
      description='A Python Code for Aggregation and Reaction',
      requires=['numpy', 'sklearn', 'scipy', 'pandas', 'matplotlib'],
      keywords='computational chemistry global minima aggregation automated reaction',
      classifiers=[
            'Development Status :: 5 - Production/Stable',
            'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
            'Programming Language :: Python :: 3.6',
            'Topic :: Scientific/Engineering :: Chemistry'
      ],
      python_requires='>=3'
      )
