"""Analysis pipeline behind the web app (no Streamlit code here, so it can also be used from scripts).

mzML -> MS2 scans -> unique precursors ("features") -> blank (background) subtraction -> adduct inference
     -> library (cosine) search + fingerprint-based candidate retrieval + formula list
"""
import os
import pickle
import sys
import threading

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import common as C  # noqa: E402

DEFAULTS = dict(
    min_peaks=5, blank_ppm=10.0, cluster_ppm=10.0, do_cluster=True,
    adduct_mode="auto",            # "auto" or a fixed adduct such as "[M+H]+"
    adduct_ppm=5.0, rt_tol=10.0,
    do_library=True, lib_ppm=20.0, frag_tol=0.02, min_matched=3, min_score=0.4, top_lib=5,
    do_retrieval=True, ret_ppm=10.0, top_ret=10,
    formula_ppm=5.0, instrument="UNKNOWN", collision_energy=None,
)


class Resources:
    """Lazy, thread-safe holder for the trained model, the spectral library and the candidate pool."""

    def __init__(self, model_path="", dataset_path=""):
        self.model_path, self.dataset_path = model_path or "", dataset_path or ""
        self._lock = threading.Lock()
        self._lib = None
        self._pool = None
        self._bundle = None
        self.library_size = 0

    @property
    def has_model(self):
        return os.path.isfile(self.model_path)

    @property
    def has_dataset(self):
        return os.path.isfile(self.dataset_path)

    def library(self):
        """(LibraryIndex, RetrievalPool) - loaded once, pool fingerprints cached next to the dataset."""
        with self._lock:
            if self._lib is None and self.has_dataset:
                with open(self.dataset_path, "rb") as f:
                    recs = pickle.load(f)
                self.library_size = len(recs)
                self._lib = C.LibraryIndex(recs)
                cache = self.dataset_path + ".pool.pkl"
                if os.path.isfile(cache):
                    self._pool = C.RetrievalPool.load(cache)
                else:
                    self._pool = C.RetrievalPool.from_records(recs)
                    try:
                        self._pool.save(cache)
                    except OSError:
                        pass
            return self._lib, self._pool

    def _model(self):
        with self._lock:
            if self._bundle is None and self.has_model:
                import torch
                import model as M
                device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
                net, ck = M.load_model(self.model_path, device)
                self._bundle = dict(
                    net=net, ck=ck, device=device, M=M, torch=torch,
                    ad={a: i for i, a in enumerate(ck["adducts"])},
                    ins={a: i for i, a in enumerate(ck.get("instruments", ["UNKNOWN"]))})
            return self._bundle

    def predict_fp(self, items):
        """items: list of (peaks, precursor_mz, adduct, instrument, ce) -> (n, 2048) probabilities."""
        b = self._model()
        if b is None:
            raise RuntimeError("model not available")
        M, torch = b["M"], b["torch"]
        out = []
        with self._lock, torch.no_grad():
            for s in range(0, len(items), 64):
                chunk = items[s:s + 64]
                x = M.make_batch([(p, mz, b["ad"].get(ad, 0), b["ins"].get(ins, 0), ce)
                                  for p, mz, ad, ins, ce in chunk], b["device"])
                _, _, fpl = b["net"].encode(*x)
                out.append(torch.sigmoid(fpl).float().cpu().numpy())
        return np.concatenate(out) if out else np.zeros((0, 2048), np.float32)


def read_candidate_list(text):
    """Parse user-supplied SMILES (one per line, extra columns ignored) into a RetrievalPool."""
    items, bad = [], 0
    for line in text.splitlines():
        tok = line.strip().split()
        if not tok or tok[0].lower() in ("smiles", "#"):
            continue
        mol, flat = C.standardize(tok[0])
        if mol is None:
            bad += 1
            continue
        info = C.mol_info(mol)
        items.append((flat, info["inchikey14"], info["exact_mass"]))
    return (C.RetrievalPool(items) if items else None), bad


def _formula_of(smiles):
    mol = C.Chem.MolFromSmiles(smiles)
    return C.rdMolDescriptors.CalcMolFormula(mol) if mol else ""


def _why_empty(name, st, min_peaks):
    if st.get("ms2", 0) == 0:
        return (f"{name}: no MS/MS (MS level 2) scans found - the file has {st.get('ms1', 0):,} MS1 scans only. "
                "Data-dependent MS/MS acquisition is needed; check the acquisition method or the msconvert settings "
                "(do not filter out MS2).")
    parts = []
    if st.get("no_precursor"):
        parts.append(f"{st['no_precursor']:,} have no precursor m/z")
    if st.get("too_few_peaks"):
        parts.append(f"{st['too_few_peaks']:,} have fewer than {min_peaks} peaks after noise removal "
                     "(lower 'Minimum peaks per spectrum' in the sidebar)")
    return f"{name}: found {st['ms2']:,} MS/MS scans but none are usable: " + "; ".join(parts) + "."


def _read(name, path, min_peaks, max_peaks, warnings=None):
    st = {}
    try:
        scans = C.read_mzml(path, min_peaks, max_peaks, stats=st)
    except C.MzMLError as exc:
        raise RuntimeError(f"{name}: {exc}") from None
    if not scans and warnings is not None:
        warnings.append(_why_empty(name, st, min_peaks))
    return scans


def run_analysis(samples, blanks, params, res, user_pool=None, progress=None):
    """samples/blanks: lists of (display_name, path_to_mzML). Returns a dict with 'features' and 'summary'."""
    P = {**DEFAULTS, **(params or {})}
    prog = progress or (lambda frac, msg: None)
    warnings = []

    # ---- blanks ----
    blank_mz = []
    for k, (name, path) in enumerate(blanks):
        prog(0.02 + 0.05 * k / max(len(blanks), 1), f"reading blank {name}")
        blank_mz += [s["precursor_mz"] for s in _read(name, path, 1, 100)]
    in_blank = C.make_blank_matcher(blank_mz, P["blank_ppm"])

    # ---- samples -> features ----
    feats, n_scans_total = [], 0
    for k, (name, path) in enumerate(samples):
        prog(0.08 + 0.12 * k / max(len(samples), 1), f"reading {name}")
        scans = _read(name, path, P["min_peaks"], 100, warnings)
        n_scans_total += len(scans)
        groups = C.cluster_scans(scans, P["cluster_ppm"]) if P["do_cluster"] else [[s] for s in scans]
        file_feats = [C.group_to_feature(g) for g in groups]
        file_feats.sort(key=lambda f: (f["rt_median"] is None, f["rt_median"] or 0.0, f["precursor_mz"]))
        for f in file_feats:
            f["file"] = name
        feats += file_feats
    for i, f in enumerate(feats, 1):
        f["id"] = f"F{i:03d}"
        f["background"] = bool(in_blank(f["precursor_mz"])) if blanks else False

    # ---- adducts (inferred per file) ----
    if P["adduct_mode"] == "auto":
        for name in {f["file"] for f in feats}:
            idx = [i for i, f in enumerate(feats) if f["file"] == name]
            inferred = C.infer_adducts([feats[i] for i in idx], P["adduct_ppm"], P["rt_tol"])
            for i, (ad, ev) in zip(idx, inferred):
                feats[i]["adduct"], feats[i]["adduct_evidence"] = ad, ev
    else:
        for f in feats:
            ok = C.ADDUCTS_BY_MODE.get(f["mode"], [])
            f["adduct"] = P["adduct_mode"] if P["adduct_mode"] in ok else ok[0]
            f["adduct_evidence"] = "fixed by user" if P["adduct_mode"] in ok else "fixed adduct not valid for this ion mode; default used"
    for f in feats:
        f["neutral_mass"] = C.neutral_mass(f["precursor_mz"], f["adduct"])
        f["formulas"] = [] if f["background"] else C.enumerate_formulas(f["neutral_mass"], P["formula_ppm"])
        f["lib_hits"], f["lib_candidates"], f["ret_hits"], f["user_hits"] = [], 0, [], []
        f["ret_candidates"] = f["user_candidates"] = 0

    work = [f for f in feats if not f["background"]]

    # ---- library search ----
    lib = pool = None
    if P["do_library"]:
        if res.has_dataset:
            prog(0.22, "loading spectral library (first run only)")
            lib, pool = res.library()
        else:
            warnings.append("No spectral library found - library search skipped.")
    if lib is not None:
        for k, f in enumerate(work):
            prog(0.35 + 0.25 * k / max(len(work), 1), f"library search {f['id']}")
            hits, ncand = lib.query(f, P["lib_ppm"], P["frag_tol"], P["min_matched"], P["top_lib"], P["min_score"])
            f["lib_candidates"] = ncand
            for rank, (s, n, r) in enumerate(hits, 1):
                f["lib_hits"].append(dict(
                    rank=rank, score=float(s), matched=int(n), smiles=r["smiles"], formula=r["formula"],
                    library_id=r["id"], adduct=r["adduct"], instrument=r["instrument"], peaks=r["peaks"],
                    confidence=C.confidence_label(s, n),
                    ppm=(r["precursor_mz"] - f["precursor_mz"]) / f["precursor_mz"] * 1e6))

    # ---- fingerprint retrieval ----
    if P["do_retrieval"] and work:
        if not res.has_model:
            warnings.append("No trained model found - fingerprint retrieval skipped.")
        else:
            prog(0.62, "predicting fingerprints")
            if pool is None and res.has_dataset:
                _, pool = res.library()
            pools = []
            if pool is not None:
                pools.append(("library", pool))
            if user_pool is not None:
                pools.append(("user_list", user_pool))
            if not pools:
                warnings.append("Retrieval needs a candidate pool (library dataset or your own SMILES list).")
            else:
                items = [(f["peaks"], f["precursor_mz"], f["adduct"], P["instrument"],
                          P["collision_energy"] if P["collision_energy"] else f.get("ce")) for f in work]
                probs = res.predict_fp(items)
                for k, (f, p) in enumerate(zip(work, probs)):
                    prog(0.75 + 0.2 * k / max(len(work), 1), f"ranking candidates {f['id']}")
                    for label, pl in pools:
                        cands, n = pl.top_candidates(p, f["neutral_mass"], P["ret_ppm"], P["top_ret"])
                        hits = [dict(rank=r, smiles=c["smiles"], formula=_formula_of(c["smiles"]), score=c["score"],
                                     ppm=(c["exact_mass"] - f["neutral_mass"]) / f["neutral_mass"] * 1e6)
                                for r, c in enumerate(cands, 1)]
                        if label == "library":
                            f["ret_hits"], f["ret_candidates"] = hits, n
                        else:
                            f["user_hits"], f["user_candidates"] = hits, n

    for f in feats:
        f["status"] = _status(f)
    prog(1.0, "done")
    return dict(features=feats, params=P, warnings=warnings, n_scans=n_scans_total,
                blank_scans=len(blank_mz), summary=summary_table(feats))


def _status(f):
    if f["background"]:
        return "Background (also in blank)"
    if f["lib_hits"]:
        c = f["lib_hits"][0]["confidence"]
        return f"Library match ({c})"
    if f["user_hits"] or f["ret_hits"]:
        return "Retrieval candidates only"
    return "No candidate (formula only)"


SUMMARY_COLS = ["ID", "File", "Precursor m/z", "RT (s)", "Scans", "Mode", "Adduct", "Adduct evidence", "Neutral mass",
                "Status", "Library SMILES", "Library cosine", "Matched peaks", "Library ID", "Retrieval top-1 SMILES",
                "Retrieval top-1 formula", "Possible formulas"]
HITS_COLS = ["ID", "File", "Precursor m/z", "RT (s)", "Adduct", "Neutral mass", "Source", "Rank", "SMILES", "Formula",
             "Score", "Matched peaks", "Mass diff (ppm)", "Library ID", "Confidence"]


def summary_table(feats):
    rows = []
    for f in feats:
        lh = f["lib_hits"][0] if f["lib_hits"] else None
        rh = (f["user_hits"] or f["ret_hits"] or [None])[0]
        rows.append({
            "ID": f["id"], "File": f["file"], "Precursor m/z": round(f["precursor_mz"], 4),
            "RT (s)": f["rt_range"], "Scans": f["n_scans"], "Mode": f["mode"],
            "Adduct": f["adduct"], "Adduct evidence": f["adduct_evidence"],
            "Neutral mass": round(f["neutral_mass"], 4), "Status": f["status"],
            "Library SMILES": lh["smiles"] if lh else "", "Library cosine": round(lh["score"], 3) if lh else None,
            "Matched peaks": lh["matched"] if lh else None, "Library ID": lh["library_id"] if lh else "",
            "Retrieval top-1 SMILES": rh["smiles"] if rh else "",
            "Retrieval top-1 formula": rh["formula"] if rh else "",
            "Possible formulas": "; ".join(f"{fm}({e:+.1f}ppm)" for fm, e in f["formulas"][:4]),
        })
    return pd.DataFrame(rows, columns=SUMMARY_COLS)


def hits_table(feats):
    rows = []
    for f in feats:
        base = {"ID": f["id"], "File": f["file"], "Precursor m/z": round(f["precursor_mz"], 5),
                "RT (s)": f["rt_range"], "Adduct": f["adduct"], "Neutral mass": round(f["neutral_mass"], 5)}
        for h in f["lib_hits"]:
            rows.append({**base, "Source": "library match", "Rank": h["rank"], "SMILES": h["smiles"],
                         "Formula": h["formula"], "Score": round(h["score"], 3), "Matched peaks": h["matched"],
                         "Mass diff (ppm)": round(h["ppm"], 1), "Library ID": h["library_id"],
                         "Confidence": h["confidence"]})
        for src, hits in (("retrieval (library molecules)", f["ret_hits"]), ("retrieval (your list)", f["user_hits"])):
            for h in hits:
                rows.append({**base, "Source": src, "Rank": h["rank"], "SMILES": h["smiles"], "Formula": h["formula"],
                             "Score": round(h["score"], 3), "Mass diff (ppm)": round(h["ppm"], 1),
                             "Confidence": "rank only"})
    return pd.DataFrame(rows, columns=HITS_COLS)
