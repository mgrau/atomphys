import pkg_resources
from pint import UnitRegistry

_ureg = UnitRegistry(system="atomic", auto_reduce_dimensions=True)
_ureg.default_format = "~0.3gP"

from .atom import Atom, elements  # noqa: E402
from .laser import Laser  # noqa: E402
from .state import State  # noqa: E402
from .transition import Transition  # noqa: E402

__version__ = pkg_resources.get_distribution("atomphys").version

__all__ = [
    "elements",
    "Atom",
    "State",
    "Transition",
    "Laser",
]
