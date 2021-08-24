def sanitize_energy(s):
    return s.strip('[]au +?')

def parse_time(f):
    if 'y' in f or 'year' in f:
        return 365.0*24.0*3600.0
    elif 'd' in f or 'day' in f:
        return 24.0*3600.0
    elif 'h' in f or 'hour' in f:
        return 3600.0
    elif 'm' in f or 'minute' in f:
        return 60.0
    elif 's' in f or 'second' in f:
        return 1.0
    elif 'ms' in f or 'millisecond' in f:
        return 1.0e-3
    elif f=='us' or f=='μs' or 'microsecond' in f:
        return 1.0e-6
    elif f=='ns' or 'nanosecond' in f:
        return 1.0e-9
    elif f=='ps' or 'picosecond' in f:
        return 1.0e-12
    elif f=='fs' or 'femtosecond' in f:
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
        x1 = x0*1.001
    fx0, fx1 = func(x0), func(x1)
    i = 2
    while (abs(fx0) > 0) and (abs((fx1-fx0)/fx0) > tol) and (i < maxfev+1):
        x2 = x1 - fx1 * (x1-x0)/(fx1 - fx0)
        x0, x1 = x1, x2
        fx0, fx1 = fx1, func(x1)
        i += 1
    return x1
