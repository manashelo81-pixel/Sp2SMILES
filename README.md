# Spectrum -> SMILES (MassBank / GNPS / MoNA)

Run in this order (all paths are examples for your machine):

```bash
pip install -r requirements.txt          # pyopenms is already installed

# 0. check how your MoNA export is laid out (reads the zip directly, no extraction)
python 00_inspect_mona.py \
  --input /home/manas/Desktop/SMILES/experimental/MoNA-export-LC-MS_Spectra-msp.zip

# 1. analyse, clean, de-duplicate, split  (GNPS .mgf and MoNA .msp / .zip supported)
python 01_prepare_data.py \
  --inputs /home/manas/Desktop/SMILES/experimental/GNPS-LIBRARY.mgf \
           /home/manas/Desktop/SMILES/experimental/MoNA-export-LC-MS_Spectra-msp.zip \
  --out /home/manas/Desktop/SMILES/processed
# GNPS and MoNA overlap heavily; duplicates are removed across both sources.
# -> dataset.pkl, dataset_meta.csv, report.json, duplicates_report.csv, dataset_overview.png

# 2. train (add --d-model 128 --layers 2 --max-train 50000 on a CPU-only laptop)
python 02_train.py --data /home/manas/Desktop/SMILES/processed/dataset.pkl \
  --out /home/manas/Desktop/SMILES/run1

# 3. predict from your own mzML files
python 03_predict.py --ckpt /home/manas/Desktop/SMILES/run1/model.pt \
  --mzml /home/manas/Desktop/SMILES/A1_GA1_01_2600.mzML --out predictions.csv
```

# 4. library search of your mzML features against the cleaned MoNA/GNPS dataset (no model needed)
python 04_library_match.py --dataset /home/manas/Desktop/SMILES/processed_all/dataset.pkl \
  --mzml /home/manas/Desktop/SMILES/A1_GA1_01_2600.mzML \
  --blank /home/manas/Desktop/SMILES/BE8-blank_BE8_01_2612.mzML /home/manas/Desktop/SMILES/BE8-blank_BE8_01_2625.mzML \
  --out library_matches.csv
```

`03_predict.py` and `04_library_match.py` both accept `--blank` (precursors also found in blanks are flagged as
background) and group repeated scans of the same precursor into one feature.

```
## What step 1 does
1. Parses spectra (streams zips, so the 3 GB MoNA files are fine).
2. Rejects: non-MS2, missing/invalid SMILES, non-organic elements, wrong size, precursor mass inconsistent
   with the structure (`--mass-tol`), too few peaks. Stereo is removed (MS rarely resolves it).
3. De-duplicates: (a) exact copies (same molecule, adduct, precursor, peaks), (b) near-duplicates within a
   molecule/adduct/mode group (cosine >= `--cos-thr`), capped at `--max-per-molecule`.
4. Splits by Murcko scaffold, so test molecules are structurally unseen.

## MoNA specifics
Beyond SMILES, the MoNA .msp parser reads precursor type (adduct), ion mode, instrument type and collision
energy. Instruments are grouped into ORBITRAP / QTOF / QQQ / ITRAP / OTHER / UNKNOWN. If a record has no SMILES
but has an InChI, the structure is rebuilt from the InChI. The model takes adduct, instrument class and collision
energy (with a "missing" flag) as inputs. The `-sdf.zip` files are not needed.

## The model
Peaks (m/z via multi-scale sinusoidal features + intensity) and the precursor/adduct feed a transformer
encoder. A decoder generates SELFIES (always valid molecules); a second head predicts a 2048-bit Morgan
fingerprint. Evaluation reports validity, exact match (InChIKey first block), mass consistency, Tanimoto,
sampled top-k with mass filtering, and fingerprint retrieval among held-out molecules.

## Caveats
- Expect low exact-match rates on unseen scaffolds; treat output as a ranked shortlist.
- Your `*-blank*.mzML` files are blanks and may have few or no MS2 spectra; the script reports how many it finds.
- mzML prediction assumes [M+H]+ (positive) or [M-H]- (negative).
- I could not read inside your MoNA zips, so the MoNA parser is written from the standard MoNA MSP layout and tested
  only on a synthetic file. Run step 0 first; if `structure_usable` is low, send me the FIRST RECORD it prints.
- Collision energy mixes eV and normalised (NCE) values; the model treats it as one number.


## Web app
`python run_webapp.py` starts a browser interface for the whole pipeline (blank subtraction, adduct inference,
library search, candidate ranking). See `webapp/README_WEB.md`.


## Deployment
See `DEPLOY.md` (GitHub + Streamlit Community Cloud) and `webapp/README_WEB.md` (local, network, Docker).
