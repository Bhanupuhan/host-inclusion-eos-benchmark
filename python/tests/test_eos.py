"""Internal consistency of the Python EoS (hostinc/eos.py)."""
import pytest
from hostinc import Diamond, Olivine, V_of, KT


@pytest.mark.parametrize('m', [Diamond(), Olivine()], ids=['diamond', 'olivine'])
def test_zero_pressure_at_reference(m):
    assert abs(m.P(m.V0, 298.0)) < 1e-12


@pytest.mark.parametrize('m', [Diamond(), Olivine()], ids=['diamond', 'olivine'])
def test_KT_equals_K0_at_reference(m):
    # at 298 K the thermal pressure is zero and (for olivine) its V-derivative too
    assert KT(m, m.V0, 298.0) == pytest.approx(m.K0, rel=1e-6)


@pytest.mark.parametrize('m', [Diamond(), Olivine()], ids=['diamond', 'olivine'])
@pytest.mark.parametrize('P,T', [(0.0, 298.0), (5.5, 1473.15), (1.0, 800.0)])
def test_volume_round_trip(m, P, T):
    assert m.P(V_of(m, P, T), T) == pytest.approx(P, abs=1e-9)
