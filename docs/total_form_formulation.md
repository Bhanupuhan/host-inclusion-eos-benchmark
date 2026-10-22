# Total-form EoS UMAT (`abaqus/umat_eos_total.f90`)

## Stress

The stress is computed from the total strain and the current temperature, never by accumulating increments:

    σ = −P_EoS(V, T) · I + 2G · dev(ε),      V = V_ref · exp(tr ε),     G = const

- ε is the total small strain from the start of the analysis (the entrapment state).
- V_ref = V(P_ref, T_ref) is found by Newton iteration on P_EoS(V, T_ref) = P_ref. It is cached per mineral.
- P_ref = 0.0055 TPa and T_ref = 1473.15 K are the last two `*User Material` constants.

At a prescribed (ε, T), this material stress evaluation has no increment-history dependence.
For the boundary-value problem, equality of final solutions also requires convergence to
the same equilibrium branch. The recorded A/B/C runs are a limited numerical check, not
a proof for arbitrary loading or geometry.

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

## Implementation scope

The reference-volume cache is keyed by mineral id and reference pressure/temperature,
not by the complete EoS parameter set. Use only one parameter set per mineral in a run.
The shared mutable cache has not been validated for threaded execution; the existing
source comment is not evidence of thread safety. Use a serial run for initial verification.

The UMAT targets the supplied 3D inputs, not a general plane-stress implementation.
The standalone compiler tests do not establish Abaqus ABI or solver compatibility.
Thermal coupling outputs are not assigned by the routine; it is intended for prescribed
temperature mechanical analyses, not coupled heat-transfer simulations.

## Literature comparison

A comparison with Moulas et al. (2023) remains for author review. No claim about that
paper's implementation is made here. The legacy increment-size test concerns only the
legacy source included in this repository.
