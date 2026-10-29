#!/usr/bin/env bash
# macOS / Linux: create a virtual environment on first use, install requirements, start the web app.
set -e
cd "$(dirname "$0")"
if [ ! -d .venv ]; then
  python3 -m venv .venv
  .venv/bin/pip install --upgrade pip
  .venv/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu
  .venv/bin/pip install -r webapp/requirements-web.txt
fi
.venv/bin/python run_webapp.py "$@"
