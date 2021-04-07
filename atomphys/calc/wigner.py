from math import log, factorial, exp
from math import floor, sqrt


def ishalfint(x):
    return 2*x == floor(2*x)


def isint(x):
    return x == floor(x)


def istriangle(a, b, c):
    return abs(a-b) <= c and c <= a + b


def Δ(a, b, c):
    return factorial(a+b-c) * factorial(a-b+c) * factorial(-a+b+c) / factorial(a+b+c+1)


wigner_3j_cache = {}
wigner_6j_cache = {}


def wigner_3j(j1, j2, j3, m1, m2, m3):
    if (j1, j2, j3, m1, m2, m3) in wigner_3j_cache:
        return wigner_3j_cache[(j1, j2, j3, m1, m2, m3)]

    if not ishalfint(j1) or not ishalfint(j2) or not ishalfint(j3) or not ishalfint(m1) or not ishalfint(m2) or not ishalfint(m3):
        raise ValueError('All arguments must be integers or half-integers')

    if (j1 > 40 or j2 > 40 or j3 > 40 or m1 > 40 or m2 > 40 or m3 > 40):
        raise OverflowError(
            'can\'t handle numbers this larger, use a sympy.physics.wigner')

    # sum of second row must equal zero
    if m1+m2+m3 != 0:
        return 0

    # triangle inequality
    if not istriangle(j1, j2, j3):
        return 0

    if j1-m1 != floor(j1-m1) or j2-m2 != floor(j2-m2) or j3-m3 != floor(j3-m3):
        return 0

    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0

    wigner = 0
    for t in range(int(max(0, j2-m1-j3, j1+m2-j3)), int(min(j1+j2-j3, j1-m1, j2+m2)) + 1):
        wigner += (-1)**t / (
            factorial(t) *
            factorial(t-(j2-m1-j3)) *
            factorial(t-(j1+m2-j3)) *
            factorial(j1+j2-j3-t) *
            factorial(j1-m1-t) *
            factorial(j2+m2-t)
        )
    wigner *= (-1)**(j1-j2-m3) * sqrt(
        Δ(j1, j2, j3) *
        factorial(j1+m1)*factorial(j1-m1) *
        factorial(j2+m2)*factorial(j2-m2) *
        factorial(j3+m3)*factorial(j3-m3)
    )

    wigner_3j_cache[(j1, j2, j3, m1, m2, m3)] = wigner

    return wigner


def wigner_6j(j1, j2, j3, J1, J2, J3):
    if (j1, j2, j3, J1, J2, J3) in wigner_6j_cache:
        return wigner_6j_cache[(j1, j2, j3, J1, J2, J3)]

    if not ishalfint(j1) or not ishalfint(j2) or not ishalfint(j3) or not ishalfint(J1) or not ishalfint(J2) or not ishalfint(J3):
        raise ValueError('All arguments must be integers or half-integers')

    # triangle inequality for each triad
    if not istriangle(j1, j2, j3) or not istriangle(j1, J2, J3) or not istriangle(J1, j2, J3) or not istriangle(J1, J2, j3):
        return 0

    # each triad must sum to an integer
    if not isint(j1+j2+j3) or not isint(j1+J2+J3) or not isint(J1+j2+J3) or not isint(J1+J2+j3):
        return 0

    wigner = 0
    for t in range(int(max(0, j1 + j2 + j3, j1 + J2 + J3, J1 + j2 + J3, J1 + J2 + j3)), int(min(j1 + j2 + J1 + J2, j2 + j3 + J2 + J3, j3 + j1 + J3 + J1))+1):
        wigner += (-1)**t * factorial(t + 1) / (
            factorial(t - (j1 + j2 + j3)) *
            factorial(t - (j1 + J2 + J3)) *
            factorial(t - (J1 + j2 + J3)) *
            factorial(t - (J1 + J2 + j3)) *
            factorial((j1 + j2 + J1 + J2) - t) *
            factorial((j2 + j3 + J2 + J3) - t) *
            factorial((j3 + j1 + J3 + J1) - t)
        )
    wigner *= sqrt(Δ(j1, j2, j3) * Δ(j1, J2, J3)
                   * Δ(J1, j2, J3) * Δ(J1, J2, j3))
    wigner_3j_cache[(j1, j2, j3, J1, J2, J3)] = wigner

    return wigner
