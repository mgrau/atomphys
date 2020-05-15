try:
    from pint import UnitRegistry
    _ureg = UnitRegistry(system='atomic')
    _Q = _ureg.Quantity
    _HAS_PINT = True
except ImportError:
    _HAS_PINT = False
    _ureg = None
    _Q = None

from .atom import *
