@echo off

REM if cl command is not found, run vcvarsall
where cl >nul 2>nul
if %ERRORLEVEL% neq 0 call vcvars32

cl /nologo /Zi /O2 /openmp main.cpp
