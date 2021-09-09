import re
from fractions import Fraction


def sanitize_energy(s):
    return s.strip('[]au +?')


def fsolve(func, x0, x1=None, tol=1.49012e-08, maxfev=100):
    '''
    Find the roots of a function.

    Return the roots of the equation ``func(x) = 0`` given a starting estimate x0

    Parameters
    ----------
    func : callable ``f(x, *args)``
        a function that takes a single argument
    x0: numeric
        The estarting estimate for the roots of ``func(x) = 0``
    tol: float, optional
        The calculation will terminate if the relative error between two consecutive iterates is at most `tol`
    maxfev: int, optional
        The maximum number of calls to the function.
    '''
    if x1 is None:
        x1 = x0 * 1.001
    fx0, fx1 = func(x0), func(x1)
    i = 2
    while (abs(fx0) > 0) and (abs((fx1 - fx0) / fx0) > tol) and (i < maxfev + 1):
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        x0, x1 = x1, x2
        fx0, fx1 = fx1, func(x1)
        i += 1
    return x1


def parse_atom_name(atom):
    atom_regex = re.search("^(\\d*)([A-Za-z]+)([0-9+]*)", atom)
    if atom_regex is not None:
        atom_isotope = atom_regex.group(1)
        atom_sym = atom_regex.group(2)
        atom_charge = atom_regex.group(3)
    else:
        raise ValueError("regex did not yield a valid atom name string")
    return atom_sym, atom_charge, atom_isotope


def parse_nuc_data(isotope):

    if 'Atomic mass' in isotope:
        if isotope['Atomic mass'] is not None:
            if len(isotope['Atomic mass']) > 0:
                try:
                    atomic_mass = float(re.search('\\d+\\.*\\d*', isotope['Atomic mass']).group())
                except AttributeError:
                    atomic_mass = isotope['Atomic mass']
            else:
                atomic_mass = 0.0
        else:
            atomic_mass = 0.0
    else:
        atomic_mass = None

    if 'Half-life' in isotope:
        if 'Stable' in isotope['Half-life']:
            half_life = 'stable'
        elif len(isotope['Half-life']) == 0:
            half_life = None
        else:
            # half_life = isotope['Half-life']
            # print(isotope['Half-life'])
            s = isotope['Half-life']
            s = re.search('^(\\d+\\.*\\d*).*\\xa0(\\w+)$', s)
            if s is not None:
                # print(s)
                value = float(s.group(1))
                t_factor = s.group(2)

                if t_factor == 'h':
                    t_factor = 'hours'
                elif t_factor == 'y':
                    t_factor = 'years'
                elif t_factor == 'μs':
                    t_factor = 'microsecond'

                half_life = value * parse_time(t_factor)
            else:
                half_life = 0.0
    else:
        half_life = None

    if 'Spin' in isotope:
        if len(isotope['Spin']) > 0:
            try:
                spin = Fraction(re.search('\\d+[/]*\\d*', isotope['Spin']).group())
            except AttributeError:
                spin = isotope['Spin (physics)']
        else:
            spin = None
    else:
        spin = None

    if 'Natural abundance' in isotope:
        if len(isotope['Natural abundance']) > 0:
            try:
                abundance = float(re.search('^\\d+\\.\\d+', isotope['Natural abundance']).group())
            except AttributeError:
                abundance = isotope['Natural abundance']
        else:
            abundance = None
    else:
        abundance = None

    if 'mag_moment_μN' in isotope:
        if isotope['mag_moment_μN'] is not None:
            try:
                value = float(re.search('\\d+\\.+\\d+', isotope['mag_moment_μN']).group())
            except AttributeError:
                value = isotope['mag_moment_μN']
            magnetic_moment = value
        else:
            magnetic_moment = 0.0
    else:
        magnetic_moment = None

    return {'spin': spin, 'gI': magnetic_moment, 'abundance': abundance, 'half_life': half_life, 'mass': atomic_mass}
