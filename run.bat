@echo off
setlocal

rem Resolve every path relative to this script's own location (not a
rem hardcoded folder), so this keeps working if/when the repo is moved or
rem re-cloned somewhere else. %~dp0 is this .bat file's directory.
set "REPO_DIR=%~dp0"
set "CODE_DIR=%REPO_DIR%code"

call "%USERPROFILE%\anaconda3\Scripts\activate.bat"
if errorlevel 1 (
    echo Could not activate Anaconda at "%USERPROFILE%\anaconda3".
    echo Edit run.bat if Anaconda is installed somewhere else.
    pause
    exit /b 1
)

rem main.py uses paths (npatlas.tsv, compoundimages\, etc.) relative to the
rem code directory, matching how Spyder's runfile() sets the working
rem directory -- so cd into code\ before launching, not the repo root.
cd /d "%CODE_DIR%"
python main.py

pause
