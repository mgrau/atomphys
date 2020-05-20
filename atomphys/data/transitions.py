import csv
import io
import urllib.request
from collections import OrderedDict

from math import pi as π


try:
    from .. import _ureg, _Q, _HAS_PINT
except ImportError:
    _HAS_PINT = False
    _ureg = None
    _Q = None


class Transition(OrderedDict):
    def __init__(self, USE_UNITS=False, **transition):
        if 'Gamma' in transition:
            Gamma = transition['Gamma']
        elif 'Aki(s^-1)' in transition:
            if USE_UNITS and _HAS_PINT:
                Gamma = _Q(
                    float(transition['Aki(s^-1)']), 's^-1').to('Eh/hbar')
            else:
                Gamma = 2.4188843265856806e-17 * float(transition['Aki(s^-1)'])
        else:
            Gamma = 0.0

        if 'Ei' in transition:
            Ei = transition['Ei']
        elif 'Ei(Ry)' in transition:
            if USE_UNITS and _HAS_PINT:
                Ei = _Q(float(transition['Ei(Ry)'].strip(
                    '[]a +?')), 'Ry').to('Eh')
            else:
                Ei = 0.5 * float(transition['Ei(Ry)'].strip('[]a +?'))
        else:
            Ei = 0.0

        if 'Ef' in transition:
            Ef = transition['Ef']
        elif 'Ek(Ry)' in transition:
            if USE_UNITS and _HAS_PINT:
                Ef = _Q(float(transition['Ek(Ry)'].strip(
                    '[]a +?')), 'Ry').to('Eh')
            else:
                Ef = 0.5 * float(transition['Ek(Ry)'].strip('[]a +?'))
        else:
            Ef = 0.0

        super(Transition, self).__init__({'Ei': Ei, 'Ef': Ef, 'Gamma': Gamma})

    def __repr__(self):
        if self.i is not None:
            state_i = '{:} {:}'.format(self.i.valence, self.i.term)
        else:
            state_i = '{:0.4g}'.format(self.Ei)
        if self.f is not None:
            state_f = '{:} {:}'.format(self.f.valence, self.f.term)
        else:
            state_f = '{:0.4g}'.format(self.Ef)

        return 'Transition({:} <---> {:}, Γ={:0.4g})'.format(state_i, state_f, self.Gamma)

    @property
    def Ei(self):
        return self['Ei']

    @property
    def Ef(self):
        return self['Ef']

    @property
    def Gamma(self):
        return self['Gamma']

    @property
    def Γ(self):
        return self['Gamma']

    @property
    def i(self):
        try:
            return self['i']
        except KeyError:
            return None

    @property
    def f(self):
        try:
            return self['f']
        except KeyError:
            return None

    @property
    def ω(self):
        if isinstance(self.Ei, _Q):
            hbar = _ureg['hbar']
        else:
            hbar = 1
        return (self.Ef-self.Ei)/hbar

    @property
    def angular_frequency(self):
        return self.ω

    @property
    def ν(self):
        return self.ω/(2*π)

    @property
    def frequency(self):
        return self.ν

    @property
    def λ(self):
        if isinstance(self.Ei, _Q):
            c = _ureg['c']
        else:
            c = 1
        return c/self.ν

    @property
    def wavelength(self):
        return self.λ


def download_nist_transitions(atom):
    # the NIST url and GET options.
    url = 'http://physics.nist.gov/cgi-bin/ASD/lines1.pl'
    values = {
        'spectra': atom,
        'format': 3,  # format {0: HTML, 1: ASCII, 2: CSV, 3: TSV}
        'en_unit': 2,  # energy units {0: cm^-1, 1: eV, 2: Ry}
        'line_out': 1,  # only with transition probabilities
        'show_av': 5,
        'allowed_out': 1,
        'forbid_out': 1,
        'enrg_out': 'on'
    }

    get_postfix = urllib.parse.urlencode(values)
    with urllib.request.urlopen(url + '?' + get_postfix) as response:
        response = response.read()

    data = csv.DictReader(io.StringIO(response.decode()), dialect='excel-tab')

    return data


def get_transitions(atom, USE_UNITS=False):
    data = download_nist_transitions(atom)
    transitions = [Transition(**row, USE_UNITS=USE_UNITS) for row in data]
    return transitions
