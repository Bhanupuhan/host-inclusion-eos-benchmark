"""hostinc: equations of state and checks for the host-inclusion EoS benchmark.

eos      - total-form EoS (BM3 + thermal pressure) for diamond and olivine, same parameters as the Abaqus UMAT
lame     - linear elastic sphere-in-sphere relaxation (independent analytical check)
results  - reader for the JOB_pinc.txt files written by abaqus/scripts/sphere_pinc.py
fem/     - PLACEHOLDER: planned Python FEM solver - not yet included
fdm/     - PLACEHOLDER: planned finite-difference solver - not yet included
"""
from .eos import EosMineral, Diamond, Olivine, V_of, KT, p_thermo
from .lame import p_inc_elastic, bulk_shear

__all__ = ['EosMineral', 'Diamond', 'Olivine', 'V_of', 'KT', 'p_thermo', 'p_inc_elastic', 'bulk_shear']
