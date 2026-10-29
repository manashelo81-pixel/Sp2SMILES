"""Analyse, clean and de-duplicate MS/MS libraries (GNPS .mgf, MoNA .msp / .zip).

Example:
  python 01_prepare_data.py \
      --inputs /home/manas/Desktop/SMILES/experimental/GNPS-LIBRARY.mgf \
      --out /home/manas/Desktop/SMILES/processed
"""
import argparse
import hashlib
import os
import zipfile
import resource
import sys
import time
import json
import pickle
import random
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import sparse

import common as C


def spectrum_vector(peaks, bin_width=0.1, n_bins=20000):
    idx = np.minimum((peaks[:, 0] / bin_width).astype(np.int64), n_bins - 1)
    v = sparse.csr_matrix((np.sqrt(peaks[:, 1]), (np.zeros(len(idx), np.int64), idx)),
                          shape=(1, n_bins))
    v.sum_duplicates()
    n = np.sqrt(v.multiply(v).sum())
    return v / n if n > 0 else v


def exact_key(r):
    p = r["peaks"]
    blob = np.round(p[:, 0], 2).tobytes() + np.round(p[:, 1] * 100).astype(np.int32).tobytes()
    return hashlib.md5(
        (r["inchikey14"] + r["adduct"] + r["mode"] + f"{r['precursor_mz']:.2f}").encode() + blob
    ).hexdigest()


def expand_inputs(paths):
    """Accept files and/or folders. Folders are searched recursively for .mgf, .msp and
    .zip files that contain .msp/.mgf. Everything else (.sdf, .crdownload, ...) is ignored."""
    found, skipped = [], []
    for p in paths:
        cands = []
        if os.path.isdir(p):
            for root, _, files in os.walk(p):
                cands += [os.path.join(root, f) for f in sorted(files)]
        else:
            cands.append(p)
        for f in cands:
            low = f.lower()
            if low.endswith((".mgf", ".msp")):
                found.append(f)
            elif low.endswith(".zip"):
                try:
                    names = zipfile.ZipFile(f).namelist()
                except Exception:
                    skipped.append((f, "unreadable zip"))
                    continue
                if any(n.lower().endswith((".msp", ".mgf")) for n in names):
                    found.append(f)
                else:
                    skipped.append((f, "zip without .msp/.mgf (e.g. SDF export)"))
            else:
                skipped.append((f, "unsupported file type"))
    return list(dict.fromkeys(found)), skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True,
                    help="files and/or folders (searched recursively for .mgf, .msp, .zip)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-peaks", type=int, default=5)
    ap.add_argument("--max-peaks", type=int, default=100)
    ap.add_argument("--max-mass", type=float, default=1000.0, help="max neutral mass (Da)")
    ap.add_argument("--mass-tol", type=float, default=0.05,
                    help="max |precursor-adduct - molecule mass| in Da (0 disables the check)")
    ap.add_argument("--cos-thr", type=float, default=0.99, help="near-duplicate cosine threshold")
    ap.add_argument("--max-per-molecule", type=int, default=10,
                    help="max kept spectra per (molecule, adduct, mode)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    import os
    os.makedirs(args.out, exist_ok=True)
    random.seed(args.seed)

    t_start = time.time()
    n_exact_removed = 0
    seen = set()
    rej = Counter()
    n_in = 0
    std_cache = {}
    kept = []

    # ---------------- 1. parse + clean ----------------
    files, skipped = expand_inputs(args.inputs)
    print("Using:")
    for f in files:
        print("  +", f)
    for f, why in skipped:
        print("  - skipped", f, f"({why})")
    if not files:
        raise SystemExit("No usable .mgf/.msp files found.")
    for path in files:
        print(f"Reading {path}", flush=True)
        for r in C.iter_records(path):
            n_in += 1
            if n_in % 50000 == 0:
                mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
                print(f"  parsed {n_in:,} spectra, kept {len(kept):,}, exact dups skipped {n_exact_removed:,}, "
                      f"{(time.time() - t_start) / 60:.1f} min, peak RAM {mem:.1f} GB", flush=True)
            if r["ms_level"] != 2:
                rej["not_ms2"] += 1
                continue
            ck = (r["smiles_raw"], r["inchi_raw"] if not r["smiles_raw"] else "")
            if ck not in std_cache:
                mol, flat = C.standardize(ck[0], inchi=ck[1])
                if mol is None:
                    std_cache[ck] = (False, flat, None, None)
                else:  # keep only light objects; RDKit Mol objects use a lot of RAM
                    std_cache[ck] = (True, flat, C.mol_info(mol), C.scaffold_of(mol))
                del mol
            ok, flat, info, scaf = std_cache[ck]
            if not ok:
                rej[flat] += 1
                continue
            if r["mode"] is None:
                rej["no_ion_mode"] += 1
                continue
            if not r["precursor_mz"]:
                rej["no_precursor"] += 1
                continue
            adduct, assumed = C.norm_adduct(r["adduct_raw"], r["name"], r["mode"])
            if adduct is None:
                rej["no_adduct"] += 1
                continue
            if info["exact_mass"] > args.max_mass:
                rej["too_heavy"] += 1
                continue
            if args.mass_tol > 0:
                err = abs(C.neutral_mass(r["precursor_mz"], adduct) - info["exact_mass"])
                if err > args.mass_tol:
                    rej["mass_mismatch"] += 1
                    continue
            peaks = C.clean_peaks(r["peaks"], r["precursor_mz"], args.max_peaks)
            if len(peaks) < args.min_peaks:
                rej["too_few_peaks"] += 1
                continue
            rec = {
                "id": r["id"], "smiles": flat, "mol_key": info["inchikey14"],
                "inchikey14": info["inchikey14"], "formula": info["formula"],
                "exact_mass": info["exact_mass"], "adduct": sys.intern(adduct), "adduct_assumed": assumed,
                "mode": sys.intern(r["mode"]), "instrument": sys.intern(r["instrument"]), "ce": r["ce"],
                "precursor_mz": float(r["precursor_mz"]), "peaks": peaks,
                "scaffold": scaf, "source": sys.intern(r["source"]),
            }
            k = exact_key(rec)           # exact duplicates are dropped immediately to save memory
            if k in seen:
                n_exact_removed += 1
                continue
            seen.add(k)
            kept.append(rec)
    print(f"Parsed {n_in:,} spectra; {len(kept):,} kept after cleaning and exact-duplicate removal")
    n_clean = len(kept) + n_exact_removed
    del std_cache, seen

    uniq = kept
    print(f"Exact duplicates removed: {n_exact_removed:,}")

    # ---------------- 3. near duplicates (cosine) + per-molecule cap ----------------
    groups = defaultdict(list)
    for i, r in enumerate(uniq):
        groups[(r["inchikey14"], r["adduct"], r["mode"])].append(i)
    keep_idx, dup_rows = [], []
    for gi, (g, idxs) in enumerate(groups.items()):
        if gi % 20000 == 0:
            print(f"  de-duplicating molecule group {gi:,}/{len(groups):,} "
                  f"({(time.time() - t_start) / 60:.1f} min)", flush=True)
        n0 = len(idxs)
        if n0 > 4000:
            idxs = random.sample(idxs, 4000)
        idxs = sorted(idxs, key=lambda i: -len(uniq[i]["peaks"]))
        if len(idxs) == 1:
            keep_idx.append(idxs[0])
        else:
            X = sparse.vstack([spectrum_vector(uniq[i]["peaks"]) for i in idxs]).tocsr()
            S = (X @ X.T).toarray()
            chosen = []
            for a in range(len(idxs)):
                if len(chosen) >= args.max_per_molecule:
                    break
                if all(S[a, b] < args.cos_thr for b in chosen):
                    chosen.append(a)
            keep_idx.extend(idxs[a] for a in chosen)
        dup_rows.append((g[0], g[1], g[2], n0, len(chosen) if len(idxs) > 1 else 1))
    final = [uniq[i] for i in sorted(keep_idx)]
    print(f"Near-duplicates / cap removed: {len(uniq) - len(final):,}; final spectra: {len(final):,}")

    # ---------------- 4. scaffold split ----------------
    for r in final:
        r["split"] = C.split_of(r["scaffold"] if r["scaffold"] else r["inchikey14"])
    split_counts = Counter(r["split"] for r in final)

    # ---------------- 5. save + analysis ----------------
    with open(f"{args.out}/dataset.pkl", "wb") as f:
        pickle.dump(final, f, protocol=4)
    df = pd.DataFrame([{k: v for k, v in r.items() if k != "peaks"} | {"n_peaks": len(r["peaks"])}
                       for r in final])
    df.to_csv(f"{args.out}/dataset_meta.csv", index=False)
    pd.DataFrame(dup_rows, columns=["inchikey14", "adduct", "mode", "n_before", "n_after"]) \
        .query("n_before != n_after").sort_values("n_before", ascending=False) \
        .to_csv(f"{args.out}/duplicates_report.csv", index=False)

    report = {
        "spectra_parsed": n_in, "rejected": dict(rej), "passed_cleaning": n_clean,
        "exact_duplicates_removed": n_exact_removed,
        "near_duplicates_or_cap_removed": len(uniq) - len(final),
        "final_spectra": len(final),
        "unique_molecules": int(df["inchikey14"].nunique()) if len(df) else 0,
        "unique_scaffolds": int(df["scaffold"].nunique()) if len(df) else 0,
        "spectra_per_split": dict(split_counts),
        "adducts": df["adduct"].value_counts().head(10).to_dict() if len(df) else {},
        "assumed_adduct_fraction": float(df["adduct_assumed"].mean()) if len(df) else 0,
        "ion_modes": df["mode"].value_counts().to_dict() if len(df) else {},
        "instrument_classes": df["instrument"].value_counts().to_dict() if len(df) else {},
        "collision_energy_available": float(df["ce"].notna().mean()) if len(df) else 0,
        "sources": df["source"].value_counts().head(10).to_dict() if len(df) else {},
        "median_peaks": float(df["n_peaks"].median()) if len(df) else 0,
    }
    with open(f"{args.out}/report.json", "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))

    if len(df):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 3, figsize=(15, 8))
        ax[0, 0].hist(df["exact_mass"], bins=60); ax[0, 0].set_title("Neutral mass (Da)")
        ax[0, 1].hist(df["n_peaks"], bins=50); ax[0, 1].set_title("Peaks per spectrum")
        df["adduct"].value_counts().head(8).plot.bar(ax=ax[0, 2]); ax[0, 2].set_title("Adducts")
        df.groupby("inchikey14").size().clip(upper=30).hist(bins=30, ax=ax[1, 0])
        ax[1, 0].set_title("Spectra per molecule (clipped at 30)")
        df["instrument"].value_counts().head(8).plot.barh(ax=ax[1, 1]); ax[1, 1].set_title("Instruments")
        pd.Series(split_counts).plot.bar(ax=ax[1, 2]); ax[1, 2].set_title("Scaffold split")
        plt.tight_layout()
        plt.savefig(f"{args.out}/dataset_overview.png", dpi=120)
    print(f"Saved to {args.out}")


if __name__ == "__main__":
    main()
