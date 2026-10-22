"""Total-form equations of state used by abaqus/umat_eos_total.f90 (same equations, same parameters).

UNITS here: pressure GPa, molar volume cm3/mol, temperature K. (The UMAT takes moduli in TPa and converts.)
Parameters: EosFit files olivine_angel2018 (Angel et al. 2018) and diamond_angel2015 (Angel et al. 2015),
reference temperature 298 K.

  isothermal part (both): Birch-Murnaghan 3rd order  P_BM3(V; V0, K0, K')
  Diamond: + Holland-Powell Einstein thermal pressure (independent of V)
      Pth = alpha0 K0 thetaE/xi0 [1/(exp(thetaE/T)-1) - 1/(exp(thetaE/298)-1)],  xi0 = u0^2 e^u0/(e^u0-1)^2, u0 = thetaE/298
  Olivine: + Mie-Grueneisen-Debye thermal pressure
      Pth = gamma(V)/V [E(T) - E(298)],  E = 3 n R T D3(theta/T),  gamma = gamma0 (V/V0)^q,
      theta = theta0 exp[(gamma0 - gamma)/q]
Diamond V0 is set to 1 (only volume ratios are used for diamond).
"""
import numpy as np
from scipy.optimize import brentq

R = 8.314462618          # J/(mol K)
TREF = 298.0             # K
_xg, _wg = np.polynomial.legendre.leggauss(40)


def D3(x):
    """Debye function D3(x) = 3/x^3 int_0^x t^3/(e^t - 1) dt (40-point Gauss-Legendre)."""
    t = x / 2 * (_xg + 1)
    return 3 / x ** 3 * np.sum(_wg * t ** 3 / np.expm1(t)) * x / 2


def bm3(V, V0, K0, Kp):
    """Birch-Murnaghan 3rd-order pressure (units of K0)."""
    f = ((V0 / V) ** (2 / 3) - 1) / 2
    return 3 * K0 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (Kp - 4) * f)


class EosMineral(object):
    """Base class: P(V, T) = P_BM3(V) + Pth(V, T). Subclasses set V0, K0, Kp and define pth(V, T)."""
    V0 = K0 = Kp = None

    def pth(self, V, T):
        raise NotImplementedError

    def P(self, V, T):
        """Pressure (GPa) at molar volume V (cm3/mol) and temperature T (K)."""
        return bm3(V, self.V0, self.K0, self.Kp) + self.pth(V, T)

    def V(self, P, T):
        """Volume at (P, T) by root finding on P(V, T) - P."""
        return V_of(self, P, T)

    def KT(self, V, T):
        """Isothermal bulk modulus K_T = -V dP/dV (GPa), central difference."""
        return KT(self, V, T)


class Diamond(EosMineral):
    K0, Kp, a0, thE = 444.0, 4.0, 2.672e-6, 1500.0     # GPa, -, 1/K, K
    V0 = 1.0

    def pth(self, V, T):
        u0 = self.thE / TREF
        xi0 = u0 * u0 * np.exp(u0) / np.expm1(u0) ** 2
        return self.a0 * self.K0 * self.thE / xi0 * (1 / np.expm1(self.thE / T) - 1 / np.expm1(self.thE / TREF))


class Olivine(EosMineral):
    K0, Kp, V0, th0, n, g0, q = 126.302, 4.54207, 43.89, 643.364, 7, 1.04424, 1.87509   # GPa, -, cm3/mol, K, -, -, -

    def pth(self, V, T):
        x = V / self.V0
        g = self.g0 * x ** self.q
        th = self.th0 * np.exp((self.g0 - g) / self.q)
        E = lambda TT: 3 * self.n * R * TT * D3(th / TT)
        return g / V * (E(T) - E(TREF)) / 1000.0          # J/cm3 = MPa -> GPa


def V_of(m, P, T):
    """Volume of mineral m at (P, T)."""
    return brentq(lambda V: m.P(V, T) - P, 0.6 * m.V0, 1.3 * m.V0, xtol=1e-14)


def KT(m, V, T, h=1e-7):
    """K_T = -V dP/dV (GPa), central difference with relative step h."""
    return -V * (m.P(V * (1 + h), T) - m.P(V * (1 - h), T)) / (2 * h * V)


def p_thermo(host, inc, P_trap, T_trap, T_end=TREF, P_end=0.0):
    """'Pthermo': inclusion pressure at (T_end) for which the inclusion's volume ratio from entrapment equals
    the host's volume ratio from (P_trap, T_trap) to (P_end, T_end), i.e. before elastic relaxation."""
    ratio = V_of(host, P_end, T_end) / V_of(host, P_trap, T_trap)
    Vi = V_of(inc, P_trap, T_trap)
    return brentq(lambda P: V_of(inc, P, T_end) / Vi - ratio, -2.0, 10.0, xtol=1e-12)
