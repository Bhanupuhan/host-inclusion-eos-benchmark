@echo off
rem Total-form EoS UMAT: build A/B/C from inputs\sphere.inp, run them, extract P_inc.  Results in work\
rem Needs: Abaqus with a Fortran compiler set up for user subroutines, and Python 3 started as "py -3".
cd /d "%~dp0"
if not exist inputs\sphere.inp (echo inputs\sphere.inp not found & exit /b 1)
if not exist umat_eos_total.f90 (echo umat_eos_total.f90 not found & exit /b 1)
if "%NCPUS%"=="" set NCPUS=4
if not exist work mkdir work
py -3 scripts\make_total.py inputs\sphere.inp work
if errorlevel 1 exit /b 1
cd work
for %%J in (sphere_total_A sphere_total_B sphere_total_C) do (
  call abaqus job=%%J input=%%J.inp user=..\umat_eos_total.f90 cpus=%NCPUS% ask_delete=OFF interactive
  call abaqus python ..\scripts\sphere_pinc.py %%J
)
type sphere_total_A_pinc.txt sphere_total_B_pinc.txt sphere_total_C_pinc.txt > TOTAL_EOS_RESULT.txt
echo written: work\TOTAL_EOS_RESULT.txt
