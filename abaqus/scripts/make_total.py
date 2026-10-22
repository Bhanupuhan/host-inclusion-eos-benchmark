"""Build the three test inputs for the TOTAL-form EoS UMAT (../umat_eos_total.f90) from the one-step input.
    py -3 make_total.py ../inputs/sphere.inp [OUT_DIR]      (source untouched)
  <stem>_total_A.inp : the original one-step setup, no trapping step
  <stem>_total_B.inp : trapping step + exhumation in 100 increments
  <stem>_total_C.inp : trapping step + exhumation in ONE increment (end state only, as EosFit)
Besides the steps, the only change is that each *User Material gets the reference state appended
(P_ref = 0.0055 TPa, T_ref = 1473.15 K): diamond 8 -> 10 constants, olivine 9 -> 11.
A path-independent material must give the same P_inc in A, B and C."""
import re
import sys
from pathlib import Path
from _common import fail, read_inp, split_step, amp_total, trap_step, write

TOOL = 'make_total.py'
P_REF, T_REF = '0.0055', '1473.15'


def props_fix(lines):
    out, i = [], 0
    while i < len(lines):
        l = lines[i]
        m = re.match(r'\*User Material,\s*constants=(\d+)\s*$', l.strip(), re.I)
        if not m:
            out.append(l)
            i += 1
            continue
        n = int(m.group(1))
        vals, comments, j = [], [], i + 1
        while len(vals) < n:
            if j >= len(lines):
                fail(TOOL, 'props of *User Material not found')
            s = lines[j].strip()
            j += 1
            if s.startswith('**'):
                comments.append(lines[j - 1])
                continue
            if not s:
                continue
            vals += [v.strip() for v in s.split(',') if v.strip()]
        if len(vals) != n:
            fail(TOOL, 'could not read %d props' % n)
        mid = int(round(float(vals[1])))
        if (mid, n) not in ((1, 8), (2, 9)):
            fail(TOOL, 'unexpected material id %d with %d props' % (mid, n))
        vals += [P_REF, T_REF]
        out.append('*User Material, constants=%d' % len(vals))
        out += comments
        out.append('** total-form EoS UMAT: last two = reference state P_ref (TPa), T_ref (K) where the strain is zero')
        for k in range(0, len(vals), 8):            # Abaqus: at most 8 values per data line
            out.append(', '.join(vals[k:k + 8]))
        i = j
    return out


if len(sys.argv) not in (2, 3):
    raise SystemExit('usage: py -3 make_total.py INPUT.inp [OUT_DIR]')
src, L = read_inp(TOOL, sys.argv[1])
out_dir = Path(sys.argv[2]) if len(sys.argv) == 3 else src.parent
out_dir.mkdir(parents=True, exist_ok=True)
L = props_fix(L)
ia, st, head, tail = split_step(TOOL, L)
trap = trap_step(st)
exh = st.replace('name=Step-1', 'name=Exh', 1)
exh1 = re.sub(r'\*Static\n[^\n]*', '*Static\n1., 1., 1., 1.', exh, count=1)


def build(tag, note, steps, with_trap):
    h = amp_total(head, ia) if with_trap else list(head)
    write(out_dir, '%s_total_%s.inp' % (src.stem, tag), h[:1] + ['** make_total.py: ' + note] + h[1:] + steps + tail, note)


build('A', 'A: original one-step setup (no trapping step)', st.split('\n'), False)
build('B', 'B: trapping step + exhumation in 100 increments', (trap + '\n' + exh).split('\n'), True)
build('C', 'C: trapping step + exhumation in ONE increment', (trap + '\n' + exh1).split('\n'), True)
