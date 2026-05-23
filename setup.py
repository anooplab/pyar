#!/usr/bin/env python3
from setuptools import setup, find_packages

INSTALL_REQUIRES = [
    'numpy>=1.18.4',
    'autograd>=1.3',
    'ase',
    'torch',
    'torchani==2.0',
    'MDAnalysis',
    'pandas>=1.0.5',
    'scipy>=1.5.2',
    'scikit-learn>=0.23.2',
    'dscribe',
    'pyh5md',
    'h5py',
    'networkx',
    'matplotlib',
    'hdbscan',
    'DBCV @ git+https://github.com/christopherjenness/DBCV.git',
    'openbabel-wheel',
]

setup(
    name='pyar',
    version='1.1.0',
    packages=find_packages(include=[
        'pyar', 'pyar.*'
    ]),
    package_data={
        'pyar': [
            'AIMNet2/models/*.jpt',
            'mlatom/MLatomF',
            'mlatom/cs.so',
            'mlatom/ref.json',
            'mlatom/aiqm1_model/*.pt',
        ]
    },
    url='https://github.com/anooplab/pyar',
    license='GPL v3',
    author='Anoop et al',
    author_email='anoop@chem.iitkgp.ac.in',
    description='A Python Code for Aggregation and Reaction',
    install_requires=INSTALL_REQUIRES,
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
            'pyar-trial-generation=pyar.scripts.trial_generation:main',
            'pyar-clustering=pyar.scripts.clustering:main',
            'pyar-benchmark-clustering=pyar.scripts.benchmark_clustering:main',
            'pyar-benchmark-orientations=pyar.scripts.benchmark_orientations:main',
            'pyar-similarity=pyar.scripts.similarity:main',
            'pyar-descriptor=pyar.scripts.descriptor:main',
            'pyar-mlopt=pyar.interface.mlopt:main',
            'pyar-aimnet2-ase-opt=pyar.AIMNet2.calculators.aimnet2_ase_opt:main',
        ],
    },
)
