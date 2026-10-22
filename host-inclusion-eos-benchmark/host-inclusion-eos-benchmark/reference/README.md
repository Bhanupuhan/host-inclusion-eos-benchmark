# Reference values

- `eosfit_reference.json`: EosFit results used as the reference (non-linear isotropic model:
  P_thermo = 0.996 GPa, P_inc = 0.842 GPa), with the system definition and the room-condition elastic
  properties printed in the same output.
- `abaqus_results.csv`: Abaqus results reported in the README.
- `raw/`: raw output backing the Abaqus values (pending).

The MATLAB script that produced the EosFit values (it calls EosFitm through MinCryst) is **not included**
until its authors agree. It is referenced in the main README. To reproduce the EosFit values, use
EosFit with the EoS files `olivine_angel2018` and `diamond_angel2015`, entrapment 5.5 GPa / 1473.15 K,
final 0 GPa / 298 K, spherical inclusion.
