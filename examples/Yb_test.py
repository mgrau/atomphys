from atomphys import _ureg, Atom
import numpy as np
import matplotlib.pyplot as plt
from math import pi as π

_ureg.default_format = '~0.3gP'
c = _ureg.c
ε_0 = _ureg.ε_0

Yb = Atom('171Yb')

# import atomphys as at
# sr_data = at.data.nuclear.get_wiki_isotope_data('Strontium','Sr')
# sr_data = at.data.nuclear.wiki_append_nuc_mag_moments(sr_data)
# sr_data[[i for i,x in enumerate(sr_data) if '87Sr' in x['Nuclide']][0]]