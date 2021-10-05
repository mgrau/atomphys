try:
    from pint import UnitRegistry
    _ureg = UnitRegistry(system='atomic', auto_reduce_dimensions=True)
    _ureg.default_format = '~0.3gP'
    _HAS_PINT = True
except ImportError:
    _HAS_PINT = False
    _ureg = None

from .atom import *
from .laser import Laser
