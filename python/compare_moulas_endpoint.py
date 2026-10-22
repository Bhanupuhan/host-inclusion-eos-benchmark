"""Evaluate Moulas et al. (2023), Eq. (3), with this repository's EoS.

Infinite spherical host, isotropic elasticity, constant host G; GPa and K.
This evaluates a published equation, not the authors' implementation.
Run from the repository root: python3 python/compare_moulas_endpoint.py
"""
import math
from scipy.optimize import brentq
from hostinc import Diamond, Olivine, V_of


def endpoint_pressures():
    host, inclusion = Diamond(), Olivine()
    p_trap, t_trap, p_final, t_final, g_host = 5.5, 1473.15, 0.0, 298.0, 531.0
    host_ratio = V_of(host, p_final, t_final) / V_of(host, p_trap, t_trap)
    inclusion_initial = V_of(inclusion, p_trap, t_trap)

    def ratio(p):
        return host_ratio * inclusion_initial / V_of(inclusion, p, t_final)

    moulas = brentq(lambda p: p - p_final - 4 * g_host * (1 - ratio(p) ** (1 / 3)), 0, 2)
    logarithmic = brentq(lambda p: p - p_final + 4 * g_host / 3 * math.log(ratio(p)), 0, 2)
    return moulas, logarithmic


if __name__ == '__main__':
    moulas, logarithmic = endpoint_pressures()
    print('Moulas 2023 Eq. (3): %.9f GPa' % moulas)
    print('Log-volume mismatch relation: %.9f GPa' % logarithmic)
    print('Reported Abaqus total-form result: 0.8409 GPa (rounded)')
    print('Relative difference, Abaqus vs Eq. (3): %.5f%%' % (100 * abs(.8409 - moulas) / moulas))
