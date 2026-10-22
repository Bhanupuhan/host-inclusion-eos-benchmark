"""Compile abaqus/umat_eos_total.f90 with gfortran and compare it with hostinc/eos.py.
Skipped when gfortran is not available. Units in the UMAT: TPa; here converted to GPa."""
import ctypes as C
import os
import shutil
import subprocess

import numpy as np
import pytest
from conftest import ROOT
from hostinc import Diamond, Olivine, V_of

DIA = [0.531, 1, 0.444, 4.0, 2.672e-6, 1500., 1.0, 0.0, 0.0055, 1473.15]
OL = [0.07583, 2, 0.126302, 4.54207, 43.89, 643.364, 7, 1.04424, 1.87509, 0.0055, 1473.15]

pytestmark = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran not found')


@pytest.fixture(scope='module')
def lib(tmp_path_factory):
    d = tmp_path_factory.mktemp('umat')
    (d / 'xit.f90').write_text("subroutine XIT\n  stop 'XIT'\nend subroutine\n")
    src = os.path.join(ROOT, 'abaqus', 'umat_eos_total.f90')
    so = str(d / 'libumat.so')
    subprocess.check_call(['gfortran', '-shared', '-fPIC', '-O2', '-J', str(d), src, str(d / 'xit.f90'), '-o', so])
    return C.CDLL(so)


def umat(lib, props, eps6, T):
    ndi, nshr = 3, 3
    nt = ndi + nshr
    d = lambda n: (C.c_double * n)()
    i = lambda v: C.byref(C.c_int(v))
    st, sv, dd = d(nt), d(3), d(nt * nt)
    z = C.c_double(0)
    lib.umat_(st, sv, dd, C.byref(z), C.byref(z), C.byref(z), C.byref(z), d(nt), d(nt), C.byref(z),
              (C.c_double * nt)(*eps6), d(nt), d(2), C.byref(z), C.byref(C.c_double(T)), C.byref(C.c_double(0)),
              d(1), d(1), C.create_string_buffer(b' ' * 80), i(ndi), i(nshr), i(nt), i(3),
              (C.c_double * len(props))(*props), i(len(props)), d(3), d(9), C.byref(z), C.byref(z), d(9), d(9),
              i(1), i(1), i(1), i(1), (C.c_int * 4)(), i(1), C.c_size_t(80))
    return np.array(st[:]), np.array(dd[:]).reshape(nt, nt).T


@pytest.mark.parametrize('props', [DIA, OL], ids=['diamond', 'olivine'])
def test_reference_state(lib, props):
    s, _ = umat(lib, props, [0.0] * 6, 1473.15)
    assert np.allclose(s[:3], -0.0055, atol=1e-12) and np.allclose(s[3:], 0.0)


@pytest.mark.parametrize('props,m', [(DIA, Diamond()), (OL, Olivine())], ids=['diamond', 'olivine'])
def test_free_crystal_at_room(lib, props, m):
    tr = np.log(V_of(m, 0.0, 298.0) / V_of(m, 5.5, 1473.15))
    s, _ = umat(lib, props, [tr / 3] * 3 + [0.0] * 3, 298.0)
    assert abs(1000 * s[0]) < 1e-6


@pytest.mark.parametrize('props,m', [(DIA, Diamond()), (OL, Olivine())], ids=['diamond', 'olivine'])
def test_constrained_cooling(lib, props, m):
    s, _ = umat(lib, props, [0.0] * 6, 298.0)
    want = m.P(V_of(m, 5.5, 1473.15), 298.0)
    assert -1000 * s[0] == pytest.approx(want, rel=1e-6)


def test_tangent_matches_finite_difference(lib):
    e0 = np.array([1e-3, -2e-3, 5e-4, 3e-4, -2e-4, 1e-4])
    s0, D = umat(lib, OL, list(e0), 800.0)
    FD = np.zeros((6, 6))
    for j in range(6):
        e = e0.copy()
        e[j] += 1e-7
        FD[:, j] = (umat(lib, OL, list(e), 800.0)[0] - s0) / 1e-7
    assert abs(D - FD).max() / abs(D).max() < 1e-4
