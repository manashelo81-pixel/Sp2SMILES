@echo off
REM Windows: create a virtual environment on first use, install requirements, start the web app.
cd /d "%~dp0"
if not exist .venv (
  python -m venv .venv
  .venv\Scripts\python -m pip install --upgrade pip
  .venv\Scripts\pip install torch --index-url https://download.pytorch.org/whl/cpu
  .venv\Scripts\pip install -r webapp\requirements-web.txt
)
.venv\Scripts\python run_webapp.py %*
