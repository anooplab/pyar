from distutils.core import setup

setup(
    name='PyAR',
    version='0.2',
    packages=['pyar', 'pyar.afir', 'pyar.data_analysis', 'pyar.interface'],
    scripts=['pyar/pyar-cli'],
    url='https://github.com/anooplab/pyar',
    license='GPl v3',
    author='anoop et al',
    author_email='anoop@chem.iitkgp.ac.in',
    description='',
    requires=['numpy', 'sklearn', 'scipy', 'pandas', 'matplotlib'],
)
