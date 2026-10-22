# host-inclusion-eos-benchmark

Verification of a non-linear, thermally consistent equation-of-state (EoS) material law for host–inclusion
elastic geobarometry in Abaqus, against EosFit.

**Benchmark.**
- An olivine sphere (r = 10 µm) in diamond.
- Entrapped at 5.5 GPa / 1473.15 K, then brought to 0 GPa / 298 K.
- Dry, no fracture.
- Reference residual inclusion pressure (EosFit, non-linear isotropic): **P_inc = 0.842 GPa**.

Authors: Biswabhanu Puhan (University of Pavia), Mattia Mazzucchelli (University of Lausanne).

## Results

All runs use the same mesh: 1/8 model, 80 000 C3D8I + 6 974 C3D10 elements.

| Method | P_inc (GPa) | vs EosFit | Notes |
|---|---|---|---|
| EosFit, non-linear isotropic (reference) | **0.842** | — | P_thermo 0.996 |
| Abaqus `*Elastic`, prestress relaxation (`sphere_prestress.inp`) | 0.8468 | +0.57 % | linear elastic route |
| Analytical elastic sphere (Lamé), same moduli and prestress | 0.8466 | +0.55 % | independent analytical check of the row above (`hostinc.lame`) |
| Abaqus, legacy incremental UMAT, no trapping step | 0.8577 | +1.9 % | superseded |
| Abaqus, legacy incremental UMAT, with trapping step | 0.7803 | −7.3 % | superseded (path dependent) |
| **Abaqus, total-form EoS UMAT** (A: no trap, B: trap + 100 increments, C: trap + 1 increment) | **0.8409** (A = B = C) | **−0.13 %** | `umat_eos_total.f90` |
| Python EoS, P_thermo only (`hostinc.eos.p_thermo`) | P_thermo 0.9957 | −0.03 % vs 0.996 | |

How to read the table:
- The Abaqus P_inc values are the **unweighted mean over the inclusion integration points**, last frame, as
  reported from the original runs.
- `abaqus/scripts/sphere_pinc.py` now also writes the IVOL volume-weighted mean. Those values need a rerun
  (the inputs now request IVOL) and are not in the table yet.
- The elastic route is compared with the Lamé solution (difference +0.03 %), not with EosFit. Its prestress of
  0.997681 GPa and its room-condition moduli make it a linear approximation of the non-linear problem.

## Repository layout

```
abaqus/
  umat_eos_total.f90        total-form EoS UMAT (main)
  legacy/                   superseded incremental UMAT + why it was replaced
  inputs/sphere.inp         3D sphere, EoS UMAT, one step (base for all EoS runs)
  inputs/sphere_prestress.inp   3D sphere, *Elastic, prestressed inclusion
  scripts/make_trap.py      adds the trapping step (legacy runs)
  scripts/make_total.py     builds the A/B/C inputs for the total-form UMAT
  scripts/sphere_pinc.py    Abaqus Python: P_inc per frame, volume-weighted and unweighted
  run_total.cmd, run_legacy.cmd, run_elastic.cmd   Windows run scripts (results in abaqus/work/)
python/
  hostinc/                  EoS (eos.py), Lamé check (lame.py), result reader (results.py)
  hostinc/fem/, hostinc/fdm/   PLACEHOLDERS for M. Mazzucchelli's Python FEM and finite-difference solver
  tests/                    pytest: EoS consistency, vs EosFit (0.5 %), recorded Abaqus results, Fortran vs Python,
                            increment-size dependence (legacy vs total form)
reference/                  EosFit reference values, Abaqus result table, raw outputs (pending)
docs/                       trapping step; total-form formulation
```

## Units

| Where | Length | Stress | Temperature | Volume |
|---|---|---|---|---|
| Abaqus inputs and UMAT | µm | **TPa** (5.5 GPa = 0.0055 TPa) | K | — |
| Python `hostinc` | — | GPa | K | cm³/mol |
| Scripts' printed output | — | GPa | K | — |

## Requirements

- Abaqus/Standard ______ (TODO: version used for the runs; the inputs were generated with Abaqus/CAE 2020),
  with user subroutines enabled.
- Fortran compiler ______ (TODO: the one used with Abaqus). The UMAT is free-form Fortran 90.
- Python 3 for the input builders (`py -3` on Windows); Abaqus Python for `sphere_pinc.py`.
- Tests: Python ≥ 3.8 with numpy, scipy and pytest (`python/requirements.txt`). gfortran is optional; without
  it the Fortran-vs-Python test is skipped.

Abaqus, EosFit and MinCryst are not included.

## How to run

```bat
cd abaqus
run_total.cmd        :: total-form UMAT, cases A, B, C  -> work\TOTAL_EOS_RESULT.txt
run_elastic.cmd      :: *Elastic prestress route        -> work\sphere_prestress_pinc.txt
run_legacy.cmd       :: legacy UMAT, no trap and trap   -> work\LEGACY_RESULT.txt
```

```bash
cd python
pip install -r requirements.txt
pytest -q tests
```

## Placeholders

- `python/hostinc/fem/` — Python FEM solver (Mattia Mazzucchelli). **Not included yet.**
- `python/hostinc/fdm/` — finite-difference solver (Mattia Mazzucchelli). **Not included yet.**
- `reference/raw/` — raw Abaqus output for each result row.

## EosFit reference

The reference values come from EosFit with:
- EoS files `olivine_angel2018` (inclusion) and `diamond_angel2015` (host);
- a spherical inclusion;
- entrapment at 5.5 GPa / 1473.15 K, final state 0 GPa / 298 K.

They were run through the authors' MATLAB script, which calls EosFitm via MinCryst (version: TODO). The script
is **not included** until its authors agree. The values and the room-condition elastic properties printed with
them are in `reference/eosfit_reference.json`.

## Changelog: changes made to the authors' files for this repository

No physics, parameters or numerical values were changed.

**`abaqus/umat_eos_total.f90`** (comments only, code unchanged):
- author line added;
- the path-independence statement now names the Abaqus checks actually run (A/B/C = 0.8409 GPa);
- EoS sources cited;
- notes added: diamond γ0 and q are read but unused; SSE/SPD/SCD are not computed; P is not clipped;
  no EoS range check; K_T is a finite difference.

**`abaqus/legacy/umat_eos_olivine_diamond.f`**:
- comments only were corrected to describe what the code does (list in `abaqus/legacy/README.md`);
- the executable code is unchanged (checked line by line).

**`abaqus/inputs/sphere.inp`** (renamed from `spheroid.inp`; the geometry is a sphere):
- olivine `*User Material` constants split over two data lines (9 values; Abaqus allows 8 per line);
- unused `*Amplitude, name=Amp-2` removed;
- unit comments corrected to TPa / µm / K;
- comments added for geometry, sets, entrapment state and SDVs;
- `IVOL` added to `*Element Output` (output only).

**`abaqus/inputs/sphere_prestress.inp`** (renamed from `spheroid2.inp`):
- unused `Amp-2` removed;
- unit and geometry comments added;
- `IVOL` added to the element output.

**Scripts** (renamed: `make_spheroid_trap.py` → `make_trap.py`, `make_spheroid_total.py` → `make_total.py`,
`spheroid_pinc.py` → `sphere_pinc.py`):
- shared code moved to `scripts/_common.py`;
- an output directory argument added; output names now follow the input name;
- `sphere_pinc.py` now also reports the IVOL volume-weighted mean.

The trap and A/B/C inputs built from the new `sphere.inp` differ from those used for the reported runs only in
the three input changes above (Amp-2, line split, IVOL). This was checked by diff.

## References

- Angel R.J., Alvaro M., Nestola F. (2018). 40 years of mineral elasticity: a critical review and a new
  parameterisation of equations of state for mantle olivines and diamond inclusions. *Physics and Chemistry of
  Minerals* 45, 95–113. https://doi.org/10.1007/s00269-017-0900-7
- Angel R.J., Alvaro M., Nestola F., Mazzucchelli M.L. (2015). Diamond thermoelastic properties and implications
  for determining the pressure of formation of diamond–inclusion systems. *Russian Geology and Geophysics* 56,
  211–220. https://doi.org/10.1016/j.rgg.2015.01.014

## Licence

Placeholder: see `LICENSE`. No licence has been granted yet.
