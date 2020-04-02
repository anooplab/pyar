from distutils.core import setup

setup(
    name='PyAR',
    version='0.2',
    packages=['pyar', 'pyar.afir', 'pyar.data_analysis', 'pyar.interface'],
    url='',
    license='GPl v3',
    author='anoop et al',
    author_email='',
    description='',
    requires=['numpy', 'sklearn', 'scipy', 'pandas', 'matplotlib'],
)
