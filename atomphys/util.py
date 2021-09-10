import re
from fractions import Fraction


def sanitize_energy(s):
    return s.strip('[]au +?')


def parse_time(f):
    if 'y' in f or 'year' in f:
        return 365.0 * 24.0 * 3600.0
    elif 'd' in f or 'day' in f:
        return 24.0 * 3600.0
    elif 'h' in f or 'hour' in f:
        return 3600.0
    elif 'm' in f or 'minute' in f:
        return 60.0
    elif 's' in f or 'second' in f:
        return 1.0
    elif 'ms' in f or 'millisecond' in f:
        return 1.0e-3
    elif f == 'us' or f == 'μs' or 'microsecond' in f:
        return 1.0e-6
    elif f == 'ns' or 'nanosecond' in f:
        return 1.0e-9
    elif f == 'ps' or 'picosecond' in f:
        return 1.0e-12
    elif f == 'fs' or 'femtosecond' in f:
        return 1.0e-15


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
    if 'Atomic mass' in isotope and isotope['Atomic mass'] is not None and len(isotope['Atomic mass']) > 0:
        try:
            atomic_mass = float(re.search('\\d+\\.*\\d*', isotope['Atomic mass']).group())
        except AttributeError:
            atomic_mass = isotope['Atomic mass']
    else:
        atomic_mass = None

    if 'Half-life' in isotope:
        if 'Stable' in isotope['Half-life']:
            half_life = 'stable'
        elif len(isotope['Half-life']) == 0:
            half_life = None
        else:
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

    if 'Spin' in isotope and len(isotope['Spin']) > 0:
        try:
            spin = float(Fraction(re.search('\\d+[/]*\\d*', isotope['Spin']).group()))
        except AttributeError:
            spin = isotope['Spin (physics)']
    else:
        spin = 0.0

    if 'Natural abundance' in isotope and len(isotope['Natural abundance']) > 0:
        try:
            abundance = float(re.search('^\\d+\\.\\d+', isotope['Natural abundance']).group())
        except AttributeError:
            abundance = isotope['Natural abundance']
    else:
        abundance = 1.0

    if 'mag_moment_μN' in isotope and isotope['mag_moment_μN'] is not None:
        try:
            value = float(re.search('\\d+\\.+\\d+', isotope['mag_moment_μN']).group())
        except AttributeError:
            value = isotope['mag_moment_μN']
        magnetic_moment = value
    else:
        magnetic_moment = 0.0

    return {'spin': spin, 'gI': magnetic_moment, 'abundance': abundance, 'half_life': half_life, 'mass': atomic_mass}

def dipole_allowed(state1,state2):
    same_selection = state1 != state2
    if 'L' in state1 and 'L' in state2:
        L_selection = abs(state1.L - state2.L) <= 1 
    else:
        L_selection = True
    J_selection = (abs(state1.J - state2.J) <= 1)
    return (same_selection and J_selection and L_selection)

def frange(start, stop=None, step=None):
    # if set start=0.0 and step = 1.0 if not specified
    start = float(start)
    if stop == None:
        stop = start + 0.0
        start = 0.0
    if step == None:
        step = 1.0
    
    num_steps = int((stop-start)/step)

    return [start + step*i for i in range(num_steps)]