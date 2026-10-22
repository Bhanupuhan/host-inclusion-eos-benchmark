# Total-form EoS UMAT (`abaqus/umat_eos_total.f90`)

## Stress

The stress is computed from the total strain and the current temperature, never by accumulating increments:

    σ = −P_EoS(V, T) · I + 2G · dev(ε),      V = V_ref · exp(tr ε),     G = const

- ε is the total small strain from the start of the analysis (the entrapment state).
- V_ref = V(P_ref, T_ref) is found by Newton iteration on P_EoS(V, T_ref) = P_ref. It is cached per mineral.
- P_ref = 0.0055 TPa and T_ref = 1473.15 K are the last two `*User Material` constants.

Because σ depends only on (ε, T), the end state cannot depend on the load path.

## EoS

The EoS is identical to `python/hostinc/eos.py`; `python/tests/test_fortran_vs_python.py` checks this.
- **Isothermal part:** Birch–Murnaghan 3rd order at 298 K.
- **Diamond:** Holland–Powell Einstein thermal pressure, independent of V.
- **Olivine:** Mie–Grüneisen–Debye thermal pressure, with γ = γ0 (V/V0)^q and
  θ = θ0 exp[(γ0 − γ)/q].
- **D3:** 20-point Gauss–Legendre quadrature in Fortran, 40-point in Python.

## Tangent

DDSDDE uses λ = K_T − 2G/3 and μ = G.
- K_T = −V dP/dV is obtained by central finite difference.
- The tangent is consistent with the stress, because dP/d(tr ε) = −K_T.

## Limitations

- G is constant.
- Small-strain kinematics (`nlgeom=NO` in the inputs); the volume ratio is taken as exp(tr ε).
- No range check on the EoS.
- SSE is not computed.

## Related literature (for the authors to confirm before publishing)

Moulas et al. (2023) state that non-linear elastic implementations for host–inclusion systems "have to be
thermodynamically consistent and have to guarantee the path-independence of the results". They also show how
incremental schemes converge with the path resolution (number of path increments). The total form used here
gives path- and increment-independent results by construction.

- Moulas E., Podladchikov Y., Zingerman K., Vershinin A., Levin V. (2023). Large-strain elastic and elasto-plastic
  formulations for host–inclusion systems and their applications in thermobarometry and geodynamics.
  *American Journal of Science* 323.
