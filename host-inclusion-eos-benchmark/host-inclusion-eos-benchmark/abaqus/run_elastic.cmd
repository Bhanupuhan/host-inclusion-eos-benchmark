@echo off
rem Linear elastic prestress route (no user subroutine). Results in work\
cd /d "%~dp0"
if not exist inputs\sphere_prestress.inp (echo inputs\sphere_prestress.inp not found & exit /b 1)
if "%NCPUS%"=="" set NCPUS=4
if not exist work mkdir work
copy /y inputs\sphere_prestress.inp work\ >nul
cd work
call abaqus job=sphere_prestress input=sphere_prestress.inp cpus=%NCPUS% ask_delete=OFF interactive
call abaqus python ..\scripts\sphere_pinc.py sphere_prestress
echo written: work\sphere_prestress_pinc.txt
