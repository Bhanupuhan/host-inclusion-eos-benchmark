"""Increment-size independence: one material point, all strain components fixed at 0 (fully constrained),
cooled from the entrapment state (5.5 GPa, 1473.15 K) to 298 K in N equal increments.
The total-form UMAT must give the same end pressure for any N, equal to the Python EoS
Needs gfortran (skipped otherwise). The ABA_PARAM.INC written here is a two-line stand-in, not the Abaqus file."""
import os
import shutil
import subprocess

import pytest
from conftest import ROOT
from hostinc import Olivine, Diamond, V_of

pytestmark = pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran not found')

DRIVER = """program inc
  implicit none
  real*8 :: st(6),sv(3),dd(6,6),ddt(6),drp(6),stran(6),dstran(6),time(2),pre(1),dpr(1),props(11)
  real*8 :: coords(3),drot(3,3),df0(3,3),df1(3,3),sse,spd,scd,rpl,drpldt,dtime,temp,dtemp,pnewdt,celent
  integer :: jstep(4),n,k,np,m
  character*80 cm
  read(*,*) m, np, n
  if (m==2) props=(/0.07583d0,2d0,0.126302d0,4.54207d0,43.89d0,643.364d0,7d0,1.04424d0,1.87509d0,0.0055d0,1473.15d0/)
  if (m==1) props=(/0.531d0,1d0,0.444d0,4d0,0.2672d-5,1500d0,1d0,0d0,0.0055d0,1473.15d0,0d0/)
  st=0d0; st(1:3)=-0.0055d0; sv=0d0; sv(1)=0.0055d0; sv(2)=1473.15d0; stran=0d0; dstran=0d0
  temp=1473.15d0; dtemp=(298d0-1473.15d0)/n
  do k=1,n
    call UMAT(st,sv,dd,sse,spd,scd,rpl,ddt,drp,drpldt,stran,dstran,time,dtime,temp,dtemp,pre,dpr,cm, &
              3,3,6,2,props,np,coords,drot,pnewdt,celent,df0,df1,1,1,1,1,jstep,1)
    temp=temp+dtemp
  end do
  write(*,*) -1000d0*sum(st(1:3))/3d0
end program
"""


@pytest.fixture(scope='module')
def exes(tmp_path_factory):
    d = tmp_path_factory.mktemp('inc')
    (d / 'drv.f90').write_text(DRIVER)
    (d / 'ABA_PARAM.INC').write_text('      IMPLICIT REAL*8(A-H,O-Z)\n      PARAMETER (NPRECD=2)\n')
    (d / 'xit.f90').write_text("subroutine XIT\n  stop 'XIT'\nend subroutine\n")
    tot = os.path.join(ROOT, 'abaqus', 'umat_eos_total.f90')
    out = {}
    for name, src in (('total', tot),):
        if not os.path.exists(src):
            pytest.skip('missing ' + src)
        exe = str(d / name)
        subprocess.check_call(['gfortran', '-O2', '-ffree-line-length-none', '-I', str(d), '-J', str(d),
                               src, str(d / 'drv.f90'), str(d / 'xit.f90'), '-o', exe], cwd=str(d))
        out[name] = exe
    return out


def run(exe, m, n):
    np_ = {('total', 1): 10, ('total', 2): 11}[(os.path.basename(exe), m)]
    return float(subprocess.run([exe], input='%d %d %d\n' % (m, np_, n), capture_output=True, text=True,
                                check=True).stdout.split()[0])


@pytest.mark.parametrize('m,mineral', [(1, Diamond()), (2, Olivine())], ids=['diamond', 'olivine'])
def test_total_form_independent_of_increments(exes, m, mineral):
    want = mineral.P(V_of(mineral, 5.5, 1473.15), 298.0)
    for n in (1, 10, 1000):
        assert run(exes['total'], m, n) == pytest.approx(want, rel=1e-6)
