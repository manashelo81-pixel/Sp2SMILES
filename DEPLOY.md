# Deploying the web app from GitHub (Streamlit Community Cloud)

## 1. What goes into GitHub, and what does not
Put the **code** in the repository (everything in this folder). Do **not** commit `model.pt` or `dataset.pkl`:
GitHub rejects files over 100 MB, and the `.gitignore` already excludes them. The app can download them at start-up
from a URL instead (step 4).

## 2. Create the repository
1. On github.com click **New repository**, name it (for example `spectrum-to-structure`), leave it empty (no README).
2. In a terminal, inside the `smiles_ms` folder:
```bash
git init
git add .
git commit -m "Spectrum to structure app"
git branch -M main
git remote add origin https://github.com/<your-user>/<your-repo>.git
git push -u origin main
```
Check `git status` first: nothing larger than 100 MB should be listed.

## 3. Deploy on Streamlit Community Cloud
1. Go to share.streamlit.io, sign in with GitHub and allow access to the repository.
2. **Create app** -> choose the repository, branch `main`, and **Main file path: `webapp/app.py`**.
3. **Advanced settings**: choose Python **3.11 or 3.12**; paste secrets (step 4).
4. **Deploy**. The first build installs `requirements.txt` and `packages.txt` (several minutes).

## 4. Model and library (optional, via secrets)
Upload `model.pt` (and optionally a library) somewhere with a direct download link, e.g. a **GitHub Release** asset
(up to 2 GB per file), Hugging Face Hub, Zenodo or cloud storage. In the app's **Settings -> Secrets**:
```toml
SMILES_MS_MODEL_URL = "https://github.com/<user>/<repo>/releases/download/v1/model.pt"
# only if the library fits in memory (see below):
# SMILES_MS_DATASET_URL = "https://github.com/<user>/<repo>/releases/download/v1/dataset.pkl"
SMILES_MS_BG = "electrons"      # electrons | spectrum | off
```
The files are downloaded once on the first visit and cached while the app runs. Without these secrets the app still
works (cleaning, blank removal, adducts, formulas; ranking from a user-supplied SMILES list needs the model).

## 5. Limits to know about
- **Memory:** Community Cloud allows about 690 MB-2.7 GB per app. The full 873,000-spectrum library needs roughly 2 GB
  or more by itself, so it will probably not fit. Use the app without the library, or a much smaller library
  (build one with `01_prepare_data.py` on a subset), or host elsewhere (below).
- Apps go to sleep when unused and wake on the next visit.
- The default upload limit is 200 MB; `.streamlit/config.toml` raises it to 500 MB (the platform may cap it).
- Public repositories are visible to everyone; the app link is public unless you restrict viewers in the settings.

## 6. If you need the full library
Use a host with more memory and run the Docker image (see `Dockerfile`): for example a Hugging Face Space (Docker),
a cloud VM, or a university server. Mount `model.pt` and `dataset.pkl` into `/app/models`.

## 7. Updating
Commit and push to `main`; Community Cloud redeploys automatically. Changing `requirements.txt` triggers a reinstall.
