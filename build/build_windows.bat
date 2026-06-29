@echo off
setlocal

rem Builds a portable, Anaconda-free Windows distribution of MPACT using
rem PyInstaller, in a throwaway venv (so it never touches the dev Anaconda
rem env). Resolves paths via %~dp0 so it works regardless of where the repo
rem is cloned. Run from anywhere; output lands in dist\MPACT\ (a folder --
rem ship the whole folder, launch via dist\MPACT\MPACT.exe).
rem
rem Requires: a plain (non-Anaconda) Python 3.9-3.12 on PATH for creating the
rem build venv. Anaconda's own Python works too if that's all you have.

set "REPO_DIR=%~dp0.."
set "VENV_DIR=%REPO_DIR%\build\.buildvenv"

if not exist "%VENV_DIR%\Scripts\python.exe" (
    echo Creating build venv at "%VENV_DIR%"...
    python -m venv "%VENV_DIR%"
    if errorlevel 1 (
        echo Failed to create venv. Is Python installed and on PATH?
        pause
        exit /b 1
    )
)

call "%VENV_DIR%\Scripts\activate.bat"

echo Installing runtime dependencies...
pip install --upgrade pip
pip install -r "%REPO_DIR%\requirements.txt"
if errorlevel 1 (
    echo Dependency install failed.
    pause
    exit /b 1
)

echo Installing PyInstaller...
pip install pyinstaller
if errorlevel 1 (
    echo PyInstaller install failed.
    pause
    exit /b 1
)

echo Building MPACT.exe...
cd /d "%REPO_DIR%"
pyinstaller build\mpact.spec --noconfirm
if errorlevel 1 (
    echo Build failed.
    pause
    exit /b 1
)

echo.
echo Build complete: %REPO_DIR%\dist\MPACT\MPACT.exe
echo Ship the entire dist\MPACT\ folder -- it is not a single-file exe.
pause
