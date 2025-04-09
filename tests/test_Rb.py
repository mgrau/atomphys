from math import pi as π

import pytest


def test_rubidium(rubidium):
    Rb = rubidium

    D2 = Rb("2S1/2").to("2P3/2")
    assert D2.ω.to("THz").m == pytest.approx(2 * π * 384.230484, rel=1e-6)
    assert D2.λ.to("nm").m == pytest.approx(780.241209, rel=1e-6)
    assert Rb("2P3/2").τ.to("ns").m == pytest.approx(26.24, rel=1e-3)
    assert D2.Γ.to("MHz").m == pytest.approx(2 * π * 6.059, rel=1e-2)

    D1 = Rb("2S1/2").to("2P1/2")
    assert D1.ω.to("THz").m == pytest.approx(2 * π * 377.107463, rel=1e-6)
    assert D1.λ.to("nm").m == pytest.approx(794.978850, rel=1e-6)
    assert Rb("2P1/2").τ.to("ns").m == pytest.approx(27.70, rel=1e-2)
    assert D1.Γ.to("MHz").m == pytest.approx(2 * π * 5.746, rel=1e-2)

    α0S = Rb("2S1/2").α()

    assert α0S.to("planck_constant Hz/(V/cm)^2").m == pytest.approx(0.0794, 4e-2)
