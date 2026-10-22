"""Shared helpers for make_total.py (input builders; the source .inp is never modified)."""
import re
from pathlib import Path


def fail(tool, m):
    raise SystemExit('%s: %s - nothing written' % (tool, m))


def read_inp(tool, path):
    src = Path(path).resolve()
    if not src.exists():
        fail(tool, 'missing ' + str(src))
    return src, src.read_text(encoding='latin-1').splitlines()


def split_step(tool, L):
    """Return (index of Amp-1 line, Step-1 text, lines before Step-1, lines after Step-1).
    Checks that the input is the expected one-step exhumation: Amp-1 = (0,1)(1,0), one *Dsload with Amp-1,
    final temperature 298 K on set united-1.all, no *Cload."""
    a = [i for i, l in enumerate(L) if re.match(r'\*Amplitude,\s*name=Amp-1\s*$', l.strip(), re.I)]
    if len(a) != 1 or [float(v) for v in L[a[0] + 1].split(',')] != [0., 1., 1., 0.]:
        fail(tool, 'Amp-1 is not the expected (0,1)(1,0) ramp')
    i0 = [i for i, l in enumerate(L) if re.match(r'\*Step,\s*name=Step-1', l, re.I)]
    if len(i0) != 1:
        fail(tool, 'Step-1 not found exactly once')
    i0 = i0[0]
    i1 = next(i for i in range(i0, len(L)) if L[i].strip().lower() == '*end step')
    if any(l.lower().startswith('*step') for l in L[i1 + 1:]):
        fail(tool, 'more than one step - unexpected')
    st = '\n'.join(L[i0:i1 + 1])
    for need in ('*Dsload, amplitude=Amp-1', 'united-1.all, 298.'):
        if st.count(need) != 1:
            fail(tool, 'Step-1 does not contain exactly one "%s"' % need)
    if '*Cload' in st:
        fail(tool, 'Step-1 has a *Cload - not handled')
    return a[0], st, L[:i0], L[i1 + 1:]


def amp_total(head, ia):
    """Amp-1 on TOTAL time: 1 during the trapping step (0-1), then 1 -> 0 during exhumation (1-2)."""
    h = list(head)
    h[ia] = '*Amplitude, name=Amp-1, time=TOTAL TIME'
    h[ia + 1] = '0., 1., 1., 1., 2., 0.'
    return h


def trap_step(st):
    """Trapping step: pressure ramped 0 -> 5.5 GPa (default ramp, no amplitude), T held at 1473.15 K.
    This balances the *Initial Conditions, type=STRESS before exhumation starts (see docs/trapping_step.md)."""
    t = st.replace('name=Step-1', 'name=Trap', 1)
    t = re.sub(r'\*Static\n[^\n]*', '*Static\n0.1, 1., 1e-05, 0.1', t, count=1)
    t = t.replace('*Dsload, amplitude=Amp-1',
                  '** trapping: load ramped 0 -> 5.5 GPa while Abaqus releases the initial-stress imbalance\n*Dsload', 1)
    return t.replace('united-1.all, 298.', 'united-1.all, 1473.15', 1)


def write(out_dir, name, lines, note):
    out = Path(out_dir) / name
    out.write_text('\r\n'.join(lines) + '\r\n', encoding='latin-1')
    print('wrote', out, '-', note)
