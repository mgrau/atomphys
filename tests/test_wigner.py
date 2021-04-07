from atomphys.calc.wigner import wigner_3j, wigner_6j
from math import sqrt
import pytest


def test_3j():
    assert wigner_3j(0, 0, 0, 0, 0, 0) == 1
    assert wigner_3j(1, 1, 0, 0, 0, 0) == pytest.approx(-sqrt(3)/3)
    assert wigner_3j(2, 1, 1, 0, 0, 0) == pytest.approx(sqrt(30)/15)
    assert wigner_3j(3, 3, 2, 2, -1, -1) == pytest.approx(sqrt(7)/14)
    assert wigner_3j(10, 6, 6, 0, 0, 0) == pytest.approx(-sqrt(73012212)/96577)
    assert wigner_3j(10, 6, 6, 10, 0, 0) == 0

    assert wigner_3j(4, 7/2, 7/2, 2, 3/2, -7/2) == pytest.approx(sqrt(165)/66)
    assert wigner_3j(7/2, 7/2, 3, 7/2, -1/2, -3) == pytest.approx(-sqrt(66)/66)
    assert wigner_3j(7/2, 7/2, 3, 7/2, 0, -3) == 0

    assert wigner_3j(4, 4, 4, 4, 0, -4) == pytest.approx(sqrt(2002)/429)
    assert wigner_3j(4, 4, 4, 4, 0, 3) == 0

    assert wigner_3j(4, 4, 9, 4, 0, -4) == 0
    assert wigner_3j(4, 9, 4, 4, 0, -4) == 0
    assert wigner_3j(9, 4, 4, 4, 0, -4) == 0

    assert wigner_3j(4, 4, 4, 5, -2, -3) == 0


def test_6j():
    assert wigner_6j(0, 0, 0, 0, 0, 0) == 1
    assert wigner_6j(2, 3/2, 3/2, 1, 3/2, 3/2) == pytest.approx(1/20)
    assert wigner_6j(3, 5/2, 1/2, 3/2, 2, 1) == pytest.approx(sqrt(30)/30)
    assert wigner_6j(7/2, 3, 3/2, 7/2, 2, 3/2) == pytest.approx(1/56)
    assert wigner_6j(4, 7/2, 3/2, 7/2, 2, 2) == pytest.approx(-sqrt(3)/84)
