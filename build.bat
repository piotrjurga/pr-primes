@echo off

set VCVARSPATH="C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars32.bat"

REM if cl command is not found, run vcvarsall
where cl >nul 2>nul
if %ERRORLEVEL% neq 0 call %VCVARSPATH%

cl /nologo /Zi /MT /O2 /Oi /openmp main.cpp
