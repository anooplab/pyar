#!/usr/bin/env python3
from setuptools import setup, find_packages

from build_backend import (
    AUTHOR,
    AUTHOR_EMAIL,
    ENTRY_POINTS,
    HOME_PAGE,
    INSTALL_REQUIRES,
    LICENSE,
    NAME,
    REQUIRES_PYTHON,
    SUMMARY,
    VERSION,
)

setup(
    name=NAME,
    version=VERSION,
    packages=find_packages(include=[
        'pyar', 'pyar.*'
    ]),
    package_data={
        'pyar': [
            'AIMNet2/models/*.jpt',
        ]
    },
    url=HOME_PAGE,
    license=LICENSE,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=SUMMARY,
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
    python_requires=REQUIRES_PYTHON,
    entry_points={
        'console_scripts': ENTRY_POINTS,
    },
)
