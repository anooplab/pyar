"""Compatibility alias for the shared growth service."""
import sys
from pyar.growth import service
sys.modules[__name__] = service
