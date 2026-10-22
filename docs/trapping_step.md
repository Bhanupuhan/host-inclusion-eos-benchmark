# The trapping step

## Problem

The one-step input (`abaqus/inputs/sphere.inp`) starts from the entrapment state through
`*Initial Conditions, type=STRESS` (−5.5 GPa in every element). It then applies the outer pressure with an
amplitude that is already at full value at t = 0: `Amp-1` = (0, 1), (1, 0).

In our own load-path tests (a 16-element block, not included here), Abaqus/Standard removed an unbalanced
initial stress linearly over the first step:
- with initial stress and **no load**, the stress fell from −5.5 GPa to 0 over the step;
- with initial stress and a **default-ramped** constant load, nothing moved;
- with initial stress and the **full-value amplitude** load, the mean stress was about 2 × the applied load
  along the path.

The end state of the step is the same either way. Only the path differs.

## Fix

The trapping step (`abaqus/scripts/make_trap.py`, `make_total.py`) has two parts:
- **Trap** (total time 0 → 1): the outer pressure is ramped 0 → 5.5 GPa with the default ramp, while T is
  held at 1473.15 K. The ramped load balances the released initial stress, so nothing moves. In the 3D sphere
  at the end of Trap, the inclusion is at 5.5004 GPa and the host at 5.4999 GPa.
- **Exh** (total time 1 → 2): `Amp-1` on total time, (0, 1), (1, 1), (2, 0). The pressure goes 5.5 → 0 GPa
  and T goes 1473.15 → 298 K.

## Why it matters

- **Path-independent material** (`*Elastic` prestress, total-form UMAT): the end value does not change.
  The total-form UMAT gives 0.8409 GPa with and without the trapping step.
- **Path-dependent material** (legacy incremental UMAT): the end value changes, 0.8577 without vs
  0.7803 GPa with the trapping step.

Open question for reviewers: we did not find this release behaviour stated explicitly in the Abaqus
documentation of `*Initial Conditions, type=STRESS`. Pointers are welcome.
