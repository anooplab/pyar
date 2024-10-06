#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name='pyar',
    version='1.0',
    packages=find_packages(include=[
        'pyar', 'pyar.*'
    ]),
    scripts=[
        'pyar/scripts/pyar-cli',
        'pyar/scripts/pyar-react',
        'pyar/scripts/pyar-explore',
        'pyar/scripts/pyar-optimiser',
        'pyar/scripts/pyar-tabu',
        'pyar/scripts/pyar-clustering',
        'pyar/scripts/pyar-similarity',
        'pyar/scripts/pyar-descriptor',
        'pyar/interface/mlopt.py',
        'pyar/AIMNet2/calculators/aimnet2_ase_opt.py'
    ],
    package_data={
        'pyar': ['AIMNet2/models/aimnet2_wb97m-d3_0.jpt']
    },
    url='https://github.com/anooplab/pyar',
    license='GPL v3',
    author='Anoop et al',
    author_email='anoop@chem.iitkgp.ac.in',
    description='A Python Code for Aggregation and Reaction',
    install_requires=[
        'numpy',
        'scikit-learn',
        'scipy',
        'pandas',
        'matplotlib',
        'pyh5md',
        'h5py',
        'hdbscan',
        'networkx',
        'DBCV @ git+https://github.com/christopherjenness/DBCV.git',
        'dscribe'
    ],
    keywords='computational chemistry global minima aggregation automated reaction',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Chemistry'
    ],
    python_requires='>=3.6',
)

