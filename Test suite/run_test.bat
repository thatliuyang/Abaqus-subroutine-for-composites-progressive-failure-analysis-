@echo off
echo =========================================
echo  Puck Unit Test - Compile and Run
echo =========================================

:: Setup Intel oneAPI Fortran environment
call "C:\Program Files (x86)\Intel\oneAPI\compiler\2025.2\env\vars.bat" >nul 2>&1

:: Compile
echo Compiling test_puck.f90 ...
ifx test_puck.f90 -o test_puck.exe 2>&1
if errorlevel 1 (
    echo.
    echo [ERROR] Compilation failed. Check the errors above.
    pause
    exit /b 1
)

echo Compile OK.
echo.

:: Run
echo Running tests...
echo.
test_puck.exe

echo.
pause
