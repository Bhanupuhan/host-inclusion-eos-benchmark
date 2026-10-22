# Legacy incremental EoS UMAT — SUPERSEDED

> **Status: superseded by `../umat_eos_total.f90`. Kept only to document why it was replaced.**
>
> `umat_eos_olivine_diamond.f` is the authors' original file. Only its **comments** were corrected
> (header, units, unused variables); the executable code is unchanged. This was checked by comparing all
> non-comment lines, with inline `!` comments removed, against the original.

## What the legacy UMAT actually does

It is an **incremental (hypoelastic)** law: σ_{n+1} = σ_n + C(K_n, G) : (Δε − α_n ΔT/3 · I).
K and α are evaluated at the start of the increment, from the stored pressure P_n (STATEV 1) and the
temperature at the start of the increment.

| | Diamond (id 1) | Olivine (id 2) |
|---|---|---|
| Bulk modulus | K = K0 + K′·P (Murnaghan-type, linear in P) | K = K0 + K′·P |
| Volume | not used | Murnaghan: V/V0 = (1 + K′P/K0)^(−1/K′) |
| Thermal expansion (volumetric) | α = α0 · ξ(T)/ξ(298), Einstein ξ | α = γ(V)·C_v(T,V) / (K(P)·V), Debye C_v |
| K(T) | none | none |
| Thermal pressure | none | none |
| G | constant | constant |

Other properties:
- P stored in STATEV(1) is clipped to ≥ 0.
- SSE is not computed.
- Units: G, K0 and P are in TPa, as in the model; there is no unit conversion except ×1e12 inside the olivine α.

## Why it was replaced

1. **Trapping-step discrepancy.** On the 3D sphere (`../inputs/sphere.inp`) it gives P_inc = **0.8577 GPa** without
   the trapping step and **0.7803 GPa** with it (EosFit 0.842).
   These recorded results differ. Establishing the cause requires checking equivalent initial states and
   increment convergence; the two values alone do not prove intrinsic material path dependence.
2. **It is not the EosFit EoS:**
   - Murnaghan K(P) instead of BM3;
   - no thermal pressure;
   - no K(T).
3. **Known open issue (not resolved):**
   - An independent integration of the same equations (a Python ODE and a 1-D spherical FE, not included) gives
     about 0.955 GPa on the trapped path.
   - Abaqus gives 0.7803.
   - The cause of this gap was not found before the law was replaced.

The total-form UMAT computes σ = −P_EoS(V, T)·I + 2G·dev(ε) from the total strain. It gives 0.8409 GPa
for all three load paths tested (see the main README).

## Increment-size dependence (reproducible: `python/tests/test_increment_dependence.py`)

Test: one material point, fully constrained (all strain components zero), cooled from 5.5 GPa / 1473.15 K
to 298 K in N equal increments. End pressure in GPa:

| N increments | 1 | 10 | 100 | 1000 | 10000 |
|---|---|---|---|---|---|
| legacy, olivine | 0.8334 | 0.9026 | 0.9333 | 0.9368 | 0.9372 |
| legacy, diamond | −2.5303 | −0.8955 | −0.6119 | −0.5824 | −0.5794 |
| total form, olivine | 0.8640 | 0.8640 | 0.8640 | 0.8640 | 0.8640 |
| total form, diamond | −0.4472 | −0.4472 | −0.4472 | −0.4472 | −0.4472 |

There are two separate defects:
1. **Explicit update.** The legacy result depends on the increment size. In Abaqus this means it also depends
   on the automatic time incrementation.
2. **A different EoS.** Even fully converged (N = 10000), the legacy law gives 0.937 instead of 0.864 for
   olivine and −0.579 instead of −0.447 for diamond. That is the effect of Murnaghan K(P), no thermal pressure
   and no K(T).

## Comment corrections made (code unchanged)

- Header rewritten to describe what the code does; marked SUPERSEDED.
- "BM3 + Holland-Powell", "BM3 + Mie-Gruneisen-Debye" and "EXACT EoS" corrected to the laws actually used.
- "G constant (path-independence requirement)" corrected: constant G alone does not establish path independence.
- GPa → TPa in the PROPS, STATEV and G comments.
- The "Convert GPa to Abaqus units … adjust as needed" comment replaced: no conversion is done.
- The "integrates alpha(T)dT correctly" comment corrected: explicit, increment dependent.
- Unused variables named in comments: DVOL, DEPS_MECH_VOL, NE_REF, NE_T, F_BM, V_RATIO, GAMMA0 and Q_PAR (in
  DIAMOND_EOS), D3_T, D3_TREF, CV_TREF, W_VAL.
- The wrong comment "(truncated for brevity - use standard tables)" removed: the 20-point table is complete.
