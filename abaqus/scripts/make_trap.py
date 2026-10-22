"""Add the trapping step to the one-step exhumation input (legacy UMAT).
    py -3 make_trap.py ../inputs/sphere.inp [OUT_DIR]   ->  OUT_DIR/sphere_trap.inp   (source untouched)
  Trap : total time 0 -> 1, pressure ramped 0 -> 5.5 GPa, T held at 1473.15 K (nothing should move)
  Exh  : total time 1 -> 2, pressure 5.5 -> 0 GPa (Amp-1 on TOTAL time), T 1473.15 -> 298 K, same increments
Same mesh, same material, same output requests."""
import sys
from pathlib import Path
from _common import read_inp, split_step, amp_total, trap_step, write

TOOL = 'make_trap.py'
if len(sys.argv) not in (2, 3):
    raise SystemExit('usage: py -3 make_trap.py INPUT.inp [OUT_DIR]')
src, L = read_inp(TOOL, sys.argv[1])
out_dir = Path(sys.argv[2]) if len(sys.argv) == 3 else src.parent
out_dir.mkdir(parents=True, exist_ok=True)
ia, st, head, tail = split_step(TOOL, L)
exh = st.replace('name=Step-1', 'name=Exh', 1)
h = amp_total(head, ia)
lines = h[:1] + ['** make_trap.py: trapping step added (Trap: equilibrium at 5.5 GPa / 1473.15 K; Exh: Amp-1 on total time)'] \
    + h[1:] + trap_step(st).split('\n') + exh.split('\n') + tail
write(out_dir, src.stem + '_trap.inp', lines, 'Trap (0->1) + Exh (1->2)')
