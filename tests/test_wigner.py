from math import sqrt

import pytest

import atomphys
from atomphys.calc.wigner import (
    ishalfint,
    isint,
    istriangle,
    wigner_3j,
    wigner_6j,
)


def test_ishalfint():
    assert ishalfint(0) is True
    assert ishalfint(1) is True
    assert ishalfint(1.0) is True
    assert ishalfint(1.5) is True
    assert ishalfint(1 / 2) is True
    assert ishalfint(1 / 3) is False
    assert ishalfint(1.50001) is False


def test_isint():
    assert isint(0) is True
    assert isint(1) is True
    assert isint(1.0) is True
    assert isint(1.5) is False
    assert isint(1 / 2) is False
    assert ishalfint(1 / 3) is False
    assert ishalfint(1.50001) is False


def test_istriangle():
    assert istriangle(1, 2, 3) is True
    assert istriangle(1, 2, 4) is False
    assert istriangle(-1, 0, 0) is False
    assert istriangle(1, -2, 2) is False
    assert istriangle(0, 0, 1) is False


def test_3j():
    assert wigner_3j(0, 0, 0, 0, 0, 0) == 1
    assert wigner_3j(1, 1, 0, 0, 0, 0) == pytest.approx(-sqrt(3) / 3)
    assert wigner_3j(2, 1, 1, 0, 0, 0) == pytest.approx(sqrt(30) / 15)
    assert wigner_3j(3, 3, 2, 2, -1, -1) == pytest.approx(sqrt(7) / 14)
    assert wigner_3j(10, 6, 6, 0, 0, 0) == pytest.approx(-sqrt(73012212) / 96577)
    assert wigner_3j(10, 6, 6, 10, 0, 0) == 0

    assert wigner_3j(4, 7 / 2, 7 / 2, 2, 3 / 2, -7 / 2) == pytest.approx(sqrt(165) / 66)
    assert wigner_3j(7 / 2, 7 / 2, 3, 7 / 2, -1 / 2, -3) == pytest.approx(
        -sqrt(66) / 66
    )
    assert wigner_3j(7 / 2, 7 / 2, 3, 7 / 2, 0, -3) == 0

    assert wigner_3j(4, 4, 4, 4, 0, -4) == pytest.approx(sqrt(2002) / 429)
    assert wigner_3j(4, 4, 4, 4, 0, 3) == 0

    assert wigner_3j(4, 4, 9, 4, 0, -4) == 0
    assert wigner_3j(4, 9, 4, 4, 0, -4) == 0
    assert wigner_3j(9, 4, 4, 4, 0, -4) == 0

    assert wigner_3j(4, 4, 4, 5, -2, -3) == 0


def test_3j_cache():
    atomphys.calc.wigner._wigner_3j_cache = {}
    assert wigner_3j(0, 0, 0, 0, 0, 0) == 1
    assert atomphys.calc.wigner._wigner_3j_cache == {(0, 0, 0, 0, 0, 0): 1.0}
    assert wigner_3j(0, 0, 0, 0, 0, 0) == 1

    # don't try this!
    atomphys.calc.wigner._wigner_3j_cache[(0, 0, 0, 0, 0, 0)] = 2
    assert wigner_3j(0, 0, 0, 0, 0, 0) == 2


def test_3j_error():
    wigner_3j(0.5, 1, 1.5, 2, 2.5, 3) == 0
    with pytest.raises(ValueError):
        wigner_3j(0.4, 1, 1.5, 2, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_3j(0.5, 1.1, 1.5, 2, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_3j(0.5, 1, 1.6, 2, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_3j(0.5, 1, 1.5, 1.9, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_3j(0.5, 1, 1.5, 2, 2.6, 3)
    with pytest.raises(ValueError):
        wigner_3j(0.5, 1, 1.5, 2, 2.5, 3 + 1e-6)

    with pytest.raises(OverflowError):
        wigner_3j(41, 0, 0, 0, 0, 0)
    with pytest.raises(OverflowError):
        wigner_3j(0, 41, 0, 0, 0, 0)
    with pytest.raises(OverflowError):
        wigner_3j(0, 0, 41, 0, 0, 0)
    with pytest.raises(OverflowError):
        wigner_3j(0, 0, 0, 41, 0, 0)
    with pytest.raises(OverflowError):
        wigner_3j(0, 0, 0, 0, 41, 0)
    with pytest.raises(OverflowError):
        wigner_3j(0, 0, 0, 0, 0, 41)


def test_wigner_3j_zero():
    assert wigner_3j(0, 0, 0, 1, 1, 1) == 0
    assert wigner_3j(0, 0, 0, 1, -1, 1) == 0
    assert wigner_3j(1, 1, 1, 1, -1, 0) == pytest.approx(sqrt(6) / 6)

    assert wigner_3j(0, 0, 0, 1 / 2, -1 / 2, 0) == 0
    assert wigner_3j(0, 0, 0, 1 / 2, 0, -1 / 2) == 0
    assert wigner_3j(0, 0, 0, 0, 1 / 2, -1 / 2) == 0


def test_6j():
    assert wigner_6j(0, 0, 0, 0, 0, 0) == 1
    assert wigner_6j(2, 3 / 2, 3 / 2, 1, 3 / 2, 3 / 2) == pytest.approx(1 / 20)
    assert wigner_6j(3, 5 / 2, 1 / 2, 3 / 2, 2, 1) == pytest.approx(sqrt(30) / 30)
    assert wigner_6j(7 / 2, 3, 3 / 2, 7 / 2, 2, 3 / 2) == pytest.approx(1 / 56)
    assert wigner_6j(4, 7 / 2, 3 / 2, 7 / 2, 2, 2) == pytest.approx(-sqrt(3) / 84)


def test_6j_cache():
    atomphys.calc.wigner._wigner_6j_cache = {}
    assert wigner_6j(0, 0, 0, 0, 0, 0) == 1
    assert atomphys.calc.wigner._wigner_6j_cache == {(0, 0, 0, 0, 0, 0): 1.0}
    assert wigner_6j(0, 0, 0, 0, 0, 0) == 1

    # don't try this!
    atomphys.calc.wigner._wigner_6j_cache[(0, 0, 0, 0, 0, 0)] = 2
    assert wigner_6j(0, 0, 0, 0, 0, 0) == 2


def test_6j_error():
    wigner_6j(0.5, 1, 1.5, 2, 2.5, 3) == 0
    with pytest.raises(ValueError):
        wigner_6j(0.4, 1, 1.5, 2, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_6j(0.5, 1.1, 1.5, 2, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_6j(0.5, 1, 1.6, 2, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_6j(0.5, 1, 1.5, 1.9, 2.5, 3)
    with pytest.raises(ValueError):
        wigner_6j(0.5, 1, 1.5, 2, 2.6, 3)
    with pytest.raises(ValueError):
        wigner_6j(0.5, 1, 1.5, 2, 2.5, 3 + 1e-6)


def test_wigner_6j_zero():
    assert wigner_6j(1, 2, 4, 0, 0, 0) == 0
    assert wigner_6j(1, 0, 0, 2, 4, 0) == 0
    assert wigner_6j(0, 1, 0, 2, 0, 4) == 0
    assert wigner_6j(0, 0, 1, 2, 4, 0) == 0

    assert wigner_6j(1 / 2, 1 / 2, 1 / 2, 1 / 2, 1 / 2, 1 / 2) == 0
