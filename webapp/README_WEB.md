# Web app

A browser interface for the whole pipeline: upload mzML files (and blanks), get grouped precursors, background
removal, adduct inference, library matches, ranked candidate structures and formulas, with spectra and structure
drawings. Results download as CSV or Excel.

## 1. Files you need
Put these in `models/` (or point to them with `--model` / `--dataset` when starting the app; the paths are not shown in the web page):
- `model.pt` - from `02_train.py` (for example your `run2/model.pt`)
- `dataset.pkl` - from `01_prepare_data.py` (the spectral library, also used as the candidate pool)

The app still starts without them: no library -> no library search; no model -> no fingerprint retrieval.

## 2. Start it
**Windows:** double-click `install_and_run.bat` (first start installs everything, a few minutes).
**macOS / Linux:** `./install_and_run.sh`
**Already have the Python environment (like your laptop):**
```bash
cd smiles_ms
pip install streamlit openpyxl      # the only new requirements
python run_webapp.py
```
A browser tab opens at http://localhost:8501. The first run with a library builds a fingerprint cache next to
`dataset.pkl` (a few minutes, once). The library needs about 2 GB of RAM.

## 3. Let other computers use it
```bash
python run_webapp.py --host 0.0.0.0 --port 8501
```
Others on your network open `http://<your-computer-ip>:8501`. The library and model are loaded once and shared by
all users. There is **no login**: do not expose it to the internet directly. For wider access put it behind a
reverse proxy with authentication (nginx, Caddy, or your institute's gateway).

## 4. Docker (any system with Docker)
```bash
docker build -t smiles-ms .
docker run -p 8501:8501 -v /path/to/models:/app/models smiles-ms
# or:  docker compose up --build
```

## 5. Notes
- Input formats: mzML and mzXML (also `.gz` compressed). Raw vendor files must be converted first (ProteoWizard msconvert).
- Upload limit is 2000 MB per file (`--max-upload-mb`).
- Use the CPU build of PyTorch unless you have a GPU: `pip install torch --index-url https://download.pytorch.org/whl/cpu`.
- Uploaded files are processed in a temporary folder and deleted when the analysis finishes.
- Candidate lists: one SMILES per line.

- The animated background (electron spirals by default) can be changed or switched off in the sidebar ("Background animation"), or at start-up with `--background off|electrons|spectrum`. It follows the light/dark setting.
