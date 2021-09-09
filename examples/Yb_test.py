from atomphys import _ureg, Atom
import numpy as np
import matplotlib.pyplot as plt
from math import pi as π

_ureg.default_format = '~0.3gP'
c = _ureg.c
ε_0 = _ureg.ε_0

Yb = Atom('Yb')