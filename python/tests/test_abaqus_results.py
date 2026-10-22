"""Recorded Abaqus results (reference/abaqus_results.csv) vs EosFit and the analytical elastic solution.
Tolerance 0.5 %. These tests check the RECORDED numbers; they do not run Abaqus."""
import pytest
from conftest import TOL, load_abaqus, load_eosfit
from hostinc import bulk_shear, p_inc_elastic

EOSFIT = load_eosfit()['results']['non_linear_isotropic']['P_inc_GPa']
AB = load_abaqus()


@pytest.mark.parametrize('case', ['total_A', 'total_B', 'total_C'])
def test_total_form_vs_eosfit(case):
    p = float(AB[case]['P_inc_ip'])
    assert abs(p - EOSFIT) / EOSFIT < TOL, (case, p, EOSFIT)


def test_recorded_total_form_cases_agree():
    ps = [float(AB[c]['P_inc_ip']) for c in ('total_A', 'total_B', 'total_C')]
    assert max(ps) - min(ps) < 1e-4


def test_elastic_vs_lame():
    # same input values as abaqus/inputs/sphere_prestress.inp (TPa)
    K_i, _ = bulk_shear(0.1896, 0.2499)
    K_h, G_h = bulk_shear(1.1391, 0.0724)
    lame = 1000.0 * p_inc_elastic(9.97681e-04, K_i, K_h, G_h)
    p = float(AB['elastic_prestress']['P_inc_ip'])
    assert abs(p - lame) / lame < TOL, (p, lame)

