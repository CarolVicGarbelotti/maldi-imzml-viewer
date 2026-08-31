@echo off
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" maldi_viewer.py
) else (
    python maldi_viewer.py
)
if errorlevel 1 pause
