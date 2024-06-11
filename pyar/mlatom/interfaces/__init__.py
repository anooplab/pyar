'''
Interfaces to third-party software.
'''

def aiqm1(**kwargs):
    from pyar.mlatom.aiqm1 import aiqm1 as interface
    return interface(**kwargs)

def ani(**kwargs):
    from pyar.mlatom.interfaces.torchani_interface import ani_methods as interface
    return interface(**kwargs)

def mndo(**kwargs):
    from pyar.mlatom.interfaces.mndo_interface import mndo_methods as interface
    return interface(**kwargs)

def sparrow(**kwargs):
    from .sparrow_interface import sparrow_methods as interface
    return interface(**kwargs)

def xtb(**kwargs):
    from pyar.mlatom.interfaces.xtb_interface import xtb_methods as interface
    return interface(**kwargs)

def dftd4(**kwargs):
    from pyar.mlatom.interfaces.dftd4_interface import dftd4_methods as interface
    return interface(**kwargs)

def ccsdtstarcbs(**kwargs):
    from pyar.mlatom.composite_methods import ccsdtstarcbs_legacy as interface
    return interface(**kwargs)

def gaussian(**kwargs):
    from pyar.mlatom.interfaces.gaussian_interface import gaussian_methods as interface
    return interface(**kwargs)

def pyscf(**kwargs):
    from pyar.mlatom.interfaces.pyscf_interface import pyscf_methods as interface
    return interface(**kwargs)