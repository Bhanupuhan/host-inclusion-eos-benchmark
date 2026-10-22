@echo off
rem LEGACY incremental UMAT (superseded): one-step input (no trapping step) and the trapped input. Results in work\
cd /d "%~dp0"
if not exist inputs\sphere.inp (echo inputs\sphere.inp not found & exit /b 1)
if not exist legacy\umat_eos_olivine_diamond.f (echo legacy\umat_eos_olivine_diamond.f not found & exit /b 1)
if "%NCPUS%"=="" set NCPUS=4
if not exist work mkdir work
copy /y inputs\sphere.inp work\ >nul
py -3 scripts\make_trap.py inputs\sphere.inp work
if errorlevel 1 exit /b 1
cd work
for %%J in (sphere sphere_trap) do (
  call abaqus job=%%J input=%%J.inp user=..\legacy\umat_eos_olivine_diamond.f cpus=%NCPUS% ask_delete=OFF interactive
  call abaqus python ..\scripts\sphere_pinc.py %%J
)
type sphere_pinc.txt sphere_trap_pinc.txt > LEGACY_RESULT.txt
echo written: work\LEGACY_RESULT.txt
