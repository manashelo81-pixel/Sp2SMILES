"""Download the trained model / spectral library from a URL on first start (for hosted deployments where the
files are too large for the code repository)."""
import os
import shutil
import tempfile
import urllib.request

CACHE_DIR = os.path.join(tempfile.gettempdir(), "smiles_ms_files")


def download(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    req = urllib.request.Request(url, headers={"User-Agent": "smiles-ms-app"})
    with urllib.request.urlopen(req, timeout=60) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 20)
    os.replace(tmp, dest)
    return dest


def ensure(path, url):
    """Return a usable path: the local file if it exists, else the downloaded copy (None url -> unchanged)."""
    if os.path.isfile(path) or not url:
        return path, None
    dest = os.path.join(CACHE_DIR, os.path.basename(path) or "download.bin")
    if os.path.isfile(dest):
        return dest, None
    try:
        return download(url, dest), None
    except Exception as exc:                       # keep the app usable without the file
        return path, f"Could not download {os.path.basename(path)}: {exc}"
