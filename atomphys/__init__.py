try:
    from pint import UnitRegistry

    _ureg = UnitRegistry(system="atomic", auto_reduce_dimensions=True)
    _ureg.default_format = "~0.3gP"
    _HAS_PINT = True
except ImportError:
    _HAS_PINT = False
    _ureg = None

import pkg_resources

from .atom import Atom, symbols
from .laser import Laser
from .states import State
from .transitions import Transition

__version__ = pkg_resources.get_distribution("atomphys").version

__all__ = [
    "symbols",
    "Atom",
    "State",
    "Transition",
    "Laser",
]
