from atomphys.util import fsolve


def test_fsolve():
    assert fsolve(lambda x: x - 5, 0) == 5
    assert fsolve(lambda x: x - 5, 1) == 5
    assert fsolve(lambda x: (x - 5) * (x + 5), 6) == 5
    assert fsolve(lambda x: (x - 5) * (x + 5), -3) == -5
    assert fsolve(lambda x: (x - 5) * (x + 5), 4, 6) == 5
    assert fsolve(lambda x: (x - 5) * (x + 5), -4, -6) == -5
