"""Compatibility adapter for the external MLatom package."""

from importlib import import_module

from pyar.optional_dependencies import optional_dependency_error

try:
    import mlatom as _mlatom
except ImportError as exc:
    raise optional_dependency_error("mlatom", feature="MLatom compatibility adapter") from exc

data = _mlatom.data
plot = _mlatom.plot
simulations = _mlatom.simulations
stats = _mlatom.stats
xyz = _mlatom.xyz
models = import_module("mlatom.models")
constants = import_module("mlatom.constants")
stopper = import_module("mlatom.stopper")
interface_MLatomF = import_module("mlatom.interface_MLatomF")
shell_cmd = import_module("mlatom.shell_cmd")

optimize_geometry = _mlatom.optimize_geometry
irc = _mlatom.irc
freq = _mlatom.freq
thermochemistry = _mlatom.thermochemistry
generate_initial_conditions = _mlatom.generate_initial_conditions
molecule = _mlatom.data.molecule

__all__ = [
    "constants",
    "data",
    "freq",
    "generate_initial_conditions",
    "interface_MLatomF",
    "irc",
    "models",
    "molecule",
    "optimize_geometry",
    "plot",
    "shell_cmd",
    "simulations",
    "stats",
    "stopper",
    "thermochemistry",
    "xyz",
]
