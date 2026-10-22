"""Linear elastic relaxation of a prestressed sphere in a spherical host (independent analytical check).

The inclusion (radius a) starts at hydrostatic pressure P0, the host is unstressed; both relax elastically.
Compatibility of the volume change at r = a:
    (P0 - P)/K_i = P * 3/(b^3 - a^3) * [a^3/(3 K_h) + b^3/(4 G_h)]       (host shell a < r < b, free at r = b)
For b -> infinity the right-hand side is 3P/(4 G_h), giving P = P0 (1/K_i) / (1/K_i + 3/(4 G_h)).
Units: any consistent stress unit.
"""


def bulk_shear(E, nu):
    """(K, G) from Young's modulus and Poisson's ratio."""
    return E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu))


def p_inc_elastic(P0, K_i, K_h, G_h, a=1.0, b=None):
    """Residual inclusion pressure after elastic relaxation. b=None: infinite host."""
    if b is None:
        c_h = 3.0 / (4.0 * G_h)
    else:
        c_h = 3.0 / (b ** 3 - a ** 3) * (a ** 3 / (3.0 * K_h) + b ** 3 / (4.0 * G_h))
    return P0 * (1.0 / K_i) / (1.0 / K_i + c_h)
