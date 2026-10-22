"""Read JOB_pinc.txt written by abaqus/scripts/sphere_pinc.py.
Returns the last frame: dict with step, t, Pinc_vw, Pinc_ip, Phst_vw, Phst_ip (GPa; nan if missing)."""
import math


def read_pinc(path):
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.replace('|', ' ').split()
            if len(parts) != 6:
                continue
            try:
                vals = [float(v) for v in parts[1:]]
            except ValueError:
                continue
            rows.append(dict(step=parts[0], t=vals[0], Pinc_vw=vals[1], Pinc_ip=vals[2], Phst_vw=vals[3], Phst_ip=vals[4]))
    if not rows:
        raise ValueError('no data rows in %s' % path)
    return rows[-1]


def rel_err(value, reference):
    return (value - reference) / reference if not math.isnan(value) else float('nan')
