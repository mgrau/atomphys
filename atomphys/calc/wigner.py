from functools import lru_cache
from math import factorial as factorial_int
from math import floor, sqrt


def factorial(x: float) -> int:
    return factorial_int(int(x))


def ishalfint(x: float) -> bool:
    """check if value is half integer

    Arguments:
        x: is this a half integer?

    Returns:
        True if x is a half integer and False if it is not.
    """
    return 2 * x == floor(2 * x)


def isint(x: float) -> bool:
    """checks if value is an integer

    Arguments:
        x: is this an integer?

    Returns:
        True if x is an integer and False if it is not.
    """
    return x == floor(x)


def istriangle(a: float, b: float, c: float) -> bool:
    """checks if triad (a, b, c) obeys the triangle inequality

    Arguments:
        a:
        b:
        c:

    Returns:
        True if the triangle inequality is satisfied and False if it is not.
    """
    return abs(a - b) <= c and c <= a + b


def Δ(a: float, b: float, c: float) -> float:
    """Helper function for wigner symbols

    Calculates the intermediate expression

    $\\Delta = \\frac{(a+b-c)!(a-b+c)!(-a+b+c)!}{(a+b+c+1)!}$

    Arguments:
        a:
        b:
        c:

    Returns:
        $\\Delta$
    """
    return (
        factorial(a + b - c)
        * factorial(a - b + c)
        * factorial(-a + b + c)
        / factorial(a + b + c + 1)
    )


@lru_cache(maxsize=None)
def wigner_3j(
    j1: float, j2: float, j3: float, m1: float, m2: float, m3: float
) -> float:
    """Calculate the Wigner 3-j symbol numerically

    The Wigner 3-j symbol is calculated numerically using a recursion relation.
    This approximate calculation is often faster than the exact analytic
    calculation that is done by `sympy.physics.wigner`.
    Due to the finite precision of floating point numbers, this is only good
    for $j$ and $m$ up to about 40, after which it becomes too inaccurate.

    The results are cached in the private variable `_wigner_3j_cache` to
    speed subsequent calculations.

    $\\left(
    \\begin{matrix}
    j_1 & j_2 & j_3 \\\\ m_1 & m_2 & m_3
    \\end{matrix}
    \\right)$

    Arguments:
        j1:
        j2:
        j3:
        m1:
        m2:
        m3:

    Returns:
        $\\left(
        \\begin{matrix}
        j_1 & j_2 & j_3 \\\\ m_1 & m_2 & m_3
        \\end{matrix}
        \\right)$
    """
    if (
        not ishalfint(j1)
        or not ishalfint(j2)
        or not ishalfint(j3)
        or not ishalfint(m1)
        or not ishalfint(m2)
        or not ishalfint(m3)
    ):
        raise ValueError("All arguments must be integers or half-integers")

    if j1 > 40 or j2 > 40 or j3 > 40 or m1 > 40 or m2 > 40 or m3 > 40:
        raise OverflowError(
            "can't handle numbers this larger, use a sympy.physics.wigner"
        )

    # sum of second row must equal zero
    if m1 + m2 + m3 != 0:
        return 0

    # triangle inequality
    if not istriangle(j1, j2, j3):
        return 0

    if (
        j1 - m1 != floor(j1 - m1)
        or j2 - m2 != floor(j2 - m2)
        or j3 - m3 != floor(j3 - m3)
    ):
        return 0

    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0

    wigner = 0
    for t in range(
        int(max(0, j2 - m1 - j3, j1 + m2 - j3)),
        int(min(j1 + j2 - j3, j1 - m1, j2 + m2)) + 1,
    ):
        wigner += (-1) ** t / (
            factorial(t)
            * factorial(t - (j2 - m1 - j3))
            * factorial(t - (j1 + m2 - j3))
            * factorial(j1 + j2 - j3 - t)
            * factorial(j1 - m1 - t)
            * factorial(j2 + m2 - t)
        )
    wigner *= (-1) ** (j1 - j2 - m3) * sqrt(
        Δ(j1, j2, j3)
        * factorial(j1 + m1)
        * factorial(j1 - m1)
        * factorial(j2 + m2)
        * factorial(j2 - m2)
        * factorial(j3 + m3)
        * factorial(j3 - m3)
    )

    return wigner


@lru_cache(maxsize=None)
def wigner_6j(
    j1: float, j2: float, j3: float, J1: float, J2: float, J3: float
) -> float:
    """Calculate the Wigner 6-j symbol numerically

    The Wigner 6-j symbol is calculated numerically using a recursion relation.
    This approximate calculation is often faster than the exact analytic
    calculation that is done by `sympy.physics.wigner`.

    The results are cached in the private variable `_wigner_6j_cache` to
    speed subsequent calculations.

    $\\left\\{
    \\begin{matrix}
    j_1 & j_2 & j_3 \\\\ J_1 & J_2 & J_3
    \\end{matrix}
    \\right\\}$

    Arguments:
        j1:
        j2:
        j3:
        J1:
        J2:
        J3:

    Returns:
        $\\left\\{
        \\begin{matrix}
        j_1 & j_2 & j_3 \\\\ J_1 & J_2 & J_3
        \\end{matrix}
        \\right\\}$
    """
    if (
        not ishalfint(j1)
        or not ishalfint(j2)
        or not ishalfint(j3)
        or not ishalfint(J1)
        or not ishalfint(J2)
        or not ishalfint(J3)
    ):
        raise ValueError("All arguments must be integers or half-integers")

    # triangle inequality for each triad
    if (
        not istriangle(j1, j2, j3)
        or not istriangle(j1, J2, J3)
        or not istriangle(J1, j2, J3)
        or not istriangle(J1, J2, j3)
    ):
        return 0

    # each triad must sum to an integer
    if (
        not isint(j1 + j2 + j3)
        or not isint(j1 + J2 + J3)
        or not isint(J1 + j2 + J3)
        or not isint(J1 + J2 + j3)
    ):
        return 0

    wigner = 0
    for t in range(
        int(max(0, j1 + j2 + j3, j1 + J2 + J3, J1 + j2 + J3, J1 + J2 + j3)),
        int(min(j1 + j2 + J1 + J2, j2 + j3 + J2 + J3, j3 + j1 + J3 + J1)) + 1,
    ):
        wigner += (
            (-1) ** t
            * factorial(t + 1)
            / (
                factorial(t - (j1 + j2 + j3))
                * factorial(t - (j1 + J2 + J3))
                * factorial(t - (J1 + j2 + J3))
                * factorial(t - (J1 + J2 + j3))
                * factorial((j1 + j2 + J1 + J2) - t)
                * factorial((j2 + j3 + J2 + J3) - t)
                * factorial((j3 + j1 + J3 + J1) - t)
            )
        )
    wigner *= sqrt(Δ(j1, j2, j3) * Δ(j1, J2, J3) * Δ(J1, j2, J3) * Δ(J1, J2, j3))

    return wigner
