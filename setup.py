#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name='pyar',
    version='1.1.0',
    packages=find_packages(include=[
        'pyar', 'pyar.*'
    ]),
    package_data={
        'pyar': [
            'AIMNet2/models/*.jpt',
            'mlatom/aiqm1_model/*.pt'
        ]
    },
    url='https://github.com/anooplab/pyar',
    license='GPL v3',
    author='Anoop et al',
    author_email='anoop@chem.iitkgp.ac.in',
    description='A Python Code for Aggregation and Reaction',
    install_requires=[
        'numpy',
        'autograd',
        'ase',
        'torch',
        'torchani',
        'MDAnalysis',
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
    entry_points={
        'console_scripts': [
            'pyar-cli=pyar.cli:main',
            'pyar-react=pyar.scripts.react:main',
            'pyar-explore=pyar.scripts.explore:main',
            'pyar-optimiser=pyar.scripts.optimiser:main',
            'pyar-tabu=pyar.scripts.tabu:main',
            'pyar-clustering=pyar.scripts.clustering:main',
            'pyar-similarity=pyar.scripts.similarity:main',
            'pyar-descriptor=pyar.scripts.descriptor:main',
            'pyar-mlopt=pyar.interface.mlopt:main',
            'pyar-aimnet2-ase-opt=pyar.AIMNet2.calculators.aimnet2_ase_opt:main',
        ],
    },
)
