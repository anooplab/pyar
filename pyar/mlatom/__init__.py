#!/usr/bin/env python3
from pyar.mlatom import data, plot, simulations, stats, xyz
from pyar.mlatom.simulations import optimize_geometry, irc, freq, thermochemistry, generate_initial_conditions
from pyar.mlatom.data import molecule

class LazyLoad:
    def __init__(self):
        self._models = None

    @property
    def models(self):
        if self._models is None:
            import importlib
            self._models = importlib.import_module('pyar.mlatom.models')
        return self._models

lazy_loader = LazyLoad()
