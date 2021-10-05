import re

re_energy = re.compile('-?\d+\.\d*|$')


def sanitize_energy(s):
    # return re_energy.findall(s)[0]

    # this is about 3.5× faster than re.findall, but it's less flexible
    # overall this can make a several hundred ms difference when loading
    return s.strip('()[]aluxyz +?').replace('&dagger;', '')


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
