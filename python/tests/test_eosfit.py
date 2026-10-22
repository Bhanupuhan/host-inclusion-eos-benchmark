"""Python EoS vs EosFit (reference/eosfit_reference.json), tolerance 0.5 %."""
from conftest import TOL, load_eosfit
from hostinc import Diamond, Olivine, p_thermo


def test_p_thermo_vs_eosfit():
    ref = load_eosfit()
    s = ref['system']
    p = p_thermo(Diamond(), Olivine(), s['entrapment']['P_GPa'], s['entrapment']['T_K'], s['final']['T_K'])
    want = ref['results']['non_linear_isotropic']['P_thermo_GPa']
    assert abs(p - want) / want < TOL, (p, want)
