"""Cross-platform launcher for the web app (Windows, macOS, Linux).

  python run_webapp.py                                  # local use: opens http://localhost:8501
  python run_webapp.py --host 0.0.0.0 --port 8501       # make it reachable from other computers on your network
  python run_webapp.py --model D:/models/model.pt --dataset D:/models/dataset.pkl
"""
import argparse
import importlib.util
import os
import subprocess
import sys
import threading
import time
import urllib.request
import webbrowser

ROOT = os.path.dirname(os.path.abspath(__file__))
REQUIRED = {"streamlit": "streamlit", "rdkit": "rdkit", "numpy": "numpy", "pandas": "pandas", "scipy": "scipy",
            "matplotlib": "matplotlib", "selfies": "selfies", "pyopenms": "pyopenms", "openpyxl": "openpyxl"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8501)
    ap.add_argument("--model", default=os.path.join(ROOT, "models", "model.pt"))
    ap.add_argument("--dataset", default=os.path.join(ROOT, "models", "dataset.pkl"))
    ap.add_argument("--lock-paths", action="store_true", help="(no longer needed; paths are never shown in the page)")
    ap.add_argument("--max-upload-mb", type=int, default=2000)
    ap.add_argument("--background", choices=["off", "electrons", "spectrum"], default="electrons",
                    help="initial page background animation (can also be changed in the sidebar)")
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()

    args.model = os.path.abspath(os.path.expanduser(args.model))       # relative paths are taken from where you ran the command
    args.dataset = os.path.abspath(os.path.expanduser(args.dataset))
    app_file = os.path.join(ROOT, "webapp", "app.py")
    needed = [os.path.join(ROOT, "webapp", n) for n in ("app.py", "pipeline.py", "glossary.py", "background.js")] + [os.path.join(ROOT, n) for n in ("common.py", "model.py")]
    absent = [p for p in needed if not os.path.isfile(p)]
    if absent:
        sys.exit("These files are missing:\n  " + "\n  ".join(absent) +
                 "\nUnzip the latest smiles_ms.zip over your smiles_ms folder (it contains the 'webapp' folder).")
    missing = [pip for mod, pip in REQUIRED.items() if importlib.util.find_spec(mod) is None]
    if missing:
        sys.exit("Missing Python packages: " + ", ".join(missing) +
                 f"\nInstall them with:  {sys.executable} -m pip install -r webapp/requirements-web.txt")
    if importlib.util.find_spec("torch") is None:
        print("Note: PyTorch is not installed - library search works, but fingerprint retrieval needs it "
              "(see webapp/README_WEB.md).")
    for label, p in (("model", args.model), ("library dataset", args.dataset)):
        if not os.path.isfile(p):
            print(f"Note: {label} not found at {p} - set it in the app's sidebar or with --{label.split()[0]}.")

    env = dict(os.environ, SMILES_MS_MODEL=args.model, SMILES_MS_DATASET=args.dataset,
               SMILES_MS_LOCK_PATHS="1" if args.lock_paths else "", SMILES_MS_BG=args.background)
    cmd = [sys.executable, "-m", "streamlit", "run", app_file,
           "--server.address", args.host, "--server.port", str(args.port),
           "--server.maxUploadSize", str(args.max_upload_mb), "--server.headless", "true",
           "--browser.gatherUsageStats", "false"]
    proc = subprocess.Popen(cmd, cwd=ROOT, env=env)
    shown = "localhost" if args.host in ("0.0.0.0", "127.0.0.1") else args.host
    url = f"http://{shown}:{args.port}"

    def wait_and_open():
        for _ in range(120):
            try:
                urllib.request.urlopen(url + "/_stcore/health", timeout=1)
                print(f"\nApp is running at {url}   (Ctrl+C to stop)")
                if not args.no_browser:
                    webbrowser.open(url)
                return
            except Exception:
                time.sleep(0.5)
    threading.Thread(target=wait_and_open, daemon=True).start()
    try:
        sys.exit(proc.wait())
    except KeyboardInterrupt:
        proc.terminate()


if __name__ == "__main__":
    main()
