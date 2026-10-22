"""abaqus python sphere_pinc.py JOB          (read-only; e.g. sphere_total_B)
Pressure P = -(S11 + S22 + S33)/3 in GPa (model stress unit TPa x 1000) for element sets INCL (inclusion) and
HST (outer host, C3D10 elements only; the inner host shell is set ENRICHMENT), for every saved frame.
Two averages over the integration points of each set:
  P_vw : volume weighted, sum(P * IVOL) / sum(IVOL)   (needs IVOL in the field output; the inputs request it)
  P_ip : simple unweighted mean over integration points (the value reported in the original runs)
Writes JOB_pinc.txt. Works with the Python 2 and Python 3 versions of Abaqus Python."""
from __future__ import print_function
import sys
import numpy as np
from odbAccess import openOdb
from abaqusConstants import INTEGRATION_POINT

if len(sys.argv) < 2:
    raise SystemExit('usage: abaqus python sphere_pinc.py JOB')
job = sys.argv[-1]
odb = openOdb(str(job + '.odb'), readOnly=True)


def find(name):
    for inst in odb.rootAssembly.instances.values():
        for k in inst.elementSets.keys():
            if k.upper() == name:
                return inst.elementSets[k]
    for k in odb.rootAssembly.elementSets.keys():
        if k.upper() == name:
            return odb.rootAssembly.elementSets[k]
    return None


def blocks(fr, var, reg):
    """(keys, data) at integration points; key = element label * 1000 + integration point number."""
    ks, ds = [], []
    for b in fr.fieldOutputs[var].getSubset(region=reg, position=INTEGRATION_POINT).bulkDataBlocks:
        lab = np.asarray(b.elementLabels, dtype=np.int64)
        ip = np.asarray(b.integrationPoints, dtype=np.int64)
        ks.append(lab * 1000 + ip)
        ds.append(np.asarray(b.data, dtype=float).reshape(len(lab), -1))
    if not ks:
        return None, None
    return np.concatenate(ks), np.concatenate(ds)


def pressures(fr, reg):
    ks, s = blocks(fr, 'S', reg)
    if ks is None:
        return float('nan'), float('nan'), 0
    p = -1000.0 * (s[:, 0] + s[:, 1] + s[:, 2]) / 3.0          # TPa -> GPa
    p_ip = float(p.mean())
    if 'IVOL' not in fr.fieldOutputs.keys():
        return float('nan'), p_ip, len(p)
    kv, v = blocks(fr, 'IVOL', reg)
    order_s, order_v = np.argsort(ks), np.argsort(kv)
    if len(ks) != len(kv) or not np.array_equal(ks[order_s], kv[order_v]):
        raise SystemExit('sphere_pinc.py: S and IVOL integration points do not match in frame %s' % fr.frameValue)
    w = v[order_v, 0]
    return float((p[order_s] * w).sum() / w.sum()), p_ip, len(p)


sets = [(n, find(n)) for n in ('INCL', 'HST')]
for n, reg in sets:
    if reg is None:
        print('warning: element set %s not found' % n)
out = ['%s   P = -mean stress (GPa). P_vw = IVOL-weighted mean, P_ip = unweighted mean over integration points' % job,
       '%-6s %8s | %9s %9s | %9s %9s' % ('step', 't', 'Pinc_vw', 'Pinc_ip', 'Phst_vw', 'Phst_ip')]
no_ivol = False
for sn in odb.steps.keys():
    for fr in odb.steps[sn].frames:
        if 'S' not in fr.fieldOutputs.keys():
            continue
        row = []
        for n, reg in sets:
            if reg is None:
                row += [float('nan'), float('nan')]
                continue
            vw, ip, _ = pressures(fr, reg)
            no_ivol = no_ivol or vw != vw
            row += [vw, ip]
        out.append('%-6s %8.4f | %9.4f %9.4f | %9.4f %9.4f' % ((sn[:6], fr.frameValue) + tuple(row)))
odb.close()
if no_ivol:
    out.append('note: IVOL missing in some frames -> P_vw = nan there (add IVOL to *Element Output and rerun)')
open(job + '_pinc.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
