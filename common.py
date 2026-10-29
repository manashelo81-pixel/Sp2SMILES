"""Shared, torch-free helpers: parsers, chemistry standardisation, peak cleaning."""
import hashlib
import io
import re
import zipfile

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

# ----------------------------------------------------------------------------
# Adducts: mass added to the neutral molecule M to get the observed m/z (z = 1)
# ----------------------------------------------------------------------------
ADDUCT_DELTA = {
    "[M+H]+": 1.007276, "[M+Na]+": 22.989218, "[M+K]+": 38.963158,
    "[M+NH4]+": 18.033823, "[M+H-H2O]+": -17.003289,
    "[M-H]-": -1.007276, "[M+Cl]-": 34.969402, "[M+FA-H]-": 44.998203,
    "[M-H-H2O]-": -19.017841,
}
_ALIASES = {}
for _k in ADDUCT_DELTA:
    _inner, _sign = _k[1:-2], _k[-1]
    for _a in (_k, _inner, _inner + _sign):
        _ALIASES[_a.lower().replace(" ", "")] = _k
_ALIASES.update({"m+hco2": "[M+FA-H]-", "[m+hcoo]-": "[M+FA-H]-", "[m-h2o+h]+": "[M+H-H2O]+"})


def norm_adduct(raw, name, mode):
    """Return (adduct, assumed_flag)."""
    if raw:
        a = _ALIASES.get(raw.strip().lower().replace(" ", ""))
        if a:
            return a, False
    if name:
        m = re.search(r"\bM[+-][A-Za-z0-9]+(?:[+-][A-Za-z0-9]+)?\b", name)
        if m:
            a = _ALIASES.get(m.group(0).lower())
            if a:
                return a, False
    if mode == "Positive":
        return "[M+H]+", True
    if mode == "Negative":
        return "[M-H]-", True
    return None, True


def norm_instrument(*texts):
    t = " ".join(x for x in texts if x).upper()
    if not t.strip():
        return "UNKNOWN"
    if any(k in t for k in ("ORBITRAP", "EXACTIVE", "Q-EXACTIVE", "FT-", "FTMS", "LTQ-FT")):
        return "ORBITRAP"
    if "TOF" in t:
        return "QTOF"
    if any(k in t for k in ("QQQ", "TRIPLE", "QTRAP", "QQ")):
        return "QQQ"
    if any(k in t for k in ("ION TRAP", "IT-", "ITMS", "LCQ", "LTQ")):
        return "ITRAP"
    return "OTHER"


def parse_ce(s):
    """First number in a collision-energy string ('35 eV', 'NCE 30', '20-40')."""
    m = re.search(r"\d+(?:\.\d+)?", s or "")
    if not m:
        return None
    v = float(m.group(0))
    return v if 0 < v <= 500 else None


def norm_mode(s):
    s = (s or "").strip().lower()
    if s.startswith("p"):
        return "Positive"
    if s.startswith("n"):
        return "Negative"
    return None


# ----------------------------------------------------------------------------
# Parsers (stream, so the multi-GB MoNA zips never need to fit in RAM)
# ----------------------------------------------------------------------------
def _to_float(x, default=None):
    try:
        return float(str(x).split()[0])
    except Exception:
        return default


def iter_mgf(lines, source=""):
    rec = None
    for line in lines:
        line = line.strip()
        if line == "BEGIN IONS":
            rec = {"meta": {}, "peaks": []}
        elif line == "END IONS":
            if rec is not None:
                yield _mgf_to_record(rec, source)
            rec = None
        elif rec is None or not line:
            continue
        elif line[0].isdigit() or line[0] in "-.":
            p = line.replace(";", " ").split()
            if len(p) >= 2:
                try:
                    rec["peaks"].append((float(p[0]), float(p[1])))
                except ValueError:
                    pass
        elif "=" in line:
            k, v = line.split("=", 1)
            rec["meta"][k.strip().upper()] = v.strip()


def _mgf_rt_seconds(m):
    for key, mult in (("RTINSECONDS", 1.0), ("RTINMINUTES", 60.0)):
        v = _to_float(m.get(key))
        if v is not None:
            return v * mult
    return None


def _mgf_to_record(rec, source):
    m = rec["meta"]
    mode = norm_mode(m.get("IONMODE") or m.get("ION_MODE"))
    if mode is None:                       # e.g. MZmine/GNPS files only carry CHARGE=1+ or 1-
        ch = m.get("CHARGE", "")
        mode = "Negative" if "-" in ch else ("Positive" if "+" in ch else None)
    return {
        "rt": _mgf_rt_seconds(m),
        "id": m.get("SPECTRUMID") or m.get("FEATURE_ID") or m.get("TITLE") or "",
        "name": m.get("NAME", ""),
        "smiles_raw": m.get("SMILES", ""),
        "inchi_raw": m.get("INCHI", ""),
        "ce": parse_ce(m.get("COLLISION_ENERGY", "")),
        "mode": mode,
        "adduct_raw": m.get("ADDUCT") or m.get("PRECURSOR_TYPE") or "",
        "precursor_mz": _to_float(m.get("PEPMASS") or m.get("PRECURSOR_MZ")),
        "instrument": norm_instrument(m.get("SOURCE_INSTRUMENT"), m.get("INSTRUMENT")),
        "ms_level": int(_to_float(m.get("MSLEVEL"), 2)),
        "peaks": np.asarray(rec["peaks"], dtype=np.float64).reshape(-1, 2),
        "source": source,
    }


def iter_msp(lines, source=""):
    """MSP as exported by MoNA. SMILES is usually inside the Comments field
    (SMILES=...); check the first record of your own export."""
    meta, peaks, npk = {}, [], None
    def flush():
        if not meta:
            return None
        com = meta.get("comments", "") + " " + meta.get("comment", "")
        smi = meta.get("smiles", "")
        if not smi:
            mm = re.search(r'SMILES[=:]\s*([^\s"]+)', com, flags=re.I)
            smi = mm.group(1) if mm else ""
        st = meta.get("spectrum_type", "MS2").upper()
        rt_min = _to_float(meta.get("retentiontime") or meta.get("retention_time") or meta.get("rt"))
        inchi = meta.get("inchi", "")
        if not inchi:
            mm = re.search(r'(InChI=1S?/[^\s"]+)', com)
            inchi = mm.group(1) if mm else ""
        return {
            "rt": rt_min * 60.0 if rt_min is not None else None,
            "id": meta.get("db#") or meta.get("accession") or meta.get("name", ""),
            "name": meta.get("name", ""),
            "smiles_raw": smi,
            "inchi_raw": inchi,
            "ce": parse_ce(meta.get("collision_energy") or meta.get("collisionenergy")
                           or meta.get("ce") or ""),
            "mode": norm_mode(meta.get("ion_mode") or meta.get("ionmode")),
            "adduct_raw": meta.get("precursor_type") or meta.get("adduct") or "",
            "precursor_mz": _to_float(meta.get("precursor_mz") or meta.get("precursormz")),
            "instrument": norm_instrument(meta.get("instrument_type"), meta.get("instrument")),
            "ms_level": 1 if st == "MS1" else 2,
            "peaks": np.asarray(peaks, dtype=np.float64).reshape(-1, 2),
            "source": source,
        }
    for line in lines:
        line = line.strip()
        if not line:
            r = flush()
            if r is not None:
                yield r
            meta, peaks, npk = {}, [], None
            continue
        if npk is not None and (line[0].isdigit() or line[0] in "-."):
            for pair in line.split(";"):
                p = pair.split()
                if len(p) >= 2:
                    try:
                        peaks.append((float(p[0]), float(p[1])))
                    except ValueError:
                        pass
        elif ":" in line:
            k, v = line.split(":", 1)
            k = k.strip().lower()
            if k == "num peaks":
                npk = int(_to_float(v, 0))
            else:
                meta[k] = v.strip()
    r = flush()
    if r is not None:
        yield r


def iter_records(path):
    p = str(path)
    low = p.lower()
    if low.endswith(".zip"):
        with zipfile.ZipFile(p) as z:
            for n in z.namelist():
                ln = n.lower()
                if ln.endswith((".msp", ".mgf")):
                    with z.open(n) as f:
                        txt = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
                        fn = iter_msp if ln.endswith(".msp") else iter_mgf
                        yield from fn(txt, source=n)
    else:
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            fn = iter_msp if low.endswith(".msp") else iter_mgf
            yield from fn(f, source=p.split("/")[-1])


# ----------------------------------------------------------------------------
# Chemistry
# ----------------------------------------------------------------------------
ALLOWED_Z = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
_FP_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def standardize(smiles, min_heavy=3, max_heavy=100, inchi=""):
    """Largest fragment, organic elements only, stereo removed. Returns
    (mol, flat_canonical_smiles) or (None, reason)."""
    s = (smiles or "").strip()
    mol = None
    if s.upper() not in {"", "N/A", "NA", "NAN", "NONE"}:
        mol = Chem.MolFromSmiles(s)
    if mol is None and inchi and inchi.startswith("InChI="):
        try:
            mol = Chem.MolFromInchi(inchi)
        except Exception:
            mol = None
    if mol is None:
        return None, ("invalid_smiles" if s and s.upper() not in {"N/A", "NA", "NAN", "NONE"} else "no_smiles")
    frags = Chem.GetMolFrags(mol, asMols=True)
    mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    if any(a.GetAtomicNum() not in ALLOWED_Z for a in mol.GetAtoms()):
        return None, "unsupported_element"
    if not (min_heavy <= mol.GetNumHeavyAtoms() <= max_heavy):
        return None, "size_filter"
    flat = Chem.MolToSmiles(mol, isomericSmiles=False)
    mol = Chem.MolFromSmiles(flat)
    if mol is None:
        return None, "invalid_smiles"
    return mol, flat


def mol_info(mol):
    return {
        "inchikey14": Chem.MolToInchiKey(mol).split("-")[0],
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "exact_mass": Descriptors.ExactMolWt(mol),
    }


def scaffold_of(mol):
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return ""


def morgan_bits(mol):
    return _FP_GEN.GetFingerprintAsNumPy(mol).astype(np.uint8)


def morgan_fp(mol):
    return _FP_GEN.GetFingerprint(mol)


def tanimoto(m1, m2):
    return DataStructs.TanimotoSimilarity(morgan_fp(m1), morgan_fp(m2))


def neutral_mass(precursor_mz, adduct):
    return precursor_mz - ADDUCT_DELTA[adduct]


def split_of(key, train=80, val=10):
    h = int(hashlib.md5(key.encode()).hexdigest(), 16) % 100
    return "train" if h < train else ("val" if h < train + val else "test")


# ----------------------------------------------------------------------------
# Peaks
# ----------------------------------------------------------------------------
def clean_peaks(peaks, precursor_mz, max_peaks=100, rel_thr=0.005):
    """Base-peak normalise, drop noise and peaks above precursor+5, keep top N,
    sort by m/z. Returns float32 (n, 2) with intensity in (0, 1]."""
    p = np.asarray(peaks, dtype=np.float64).reshape(-1, 2)
    p = p[(p[:, 1] > 0) & (p[:, 0] > 0)]
    if precursor_mz:
        p = p[p[:, 0] <= precursor_mz + 5.0]
    if len(p) == 0:
        return np.zeros((0, 2), np.float32)
    p[:, 1] /= p[:, 1].max()
    p = p[p[:, 1] >= rel_thr]
    if len(p) > max_peaks:
        p = p[np.argsort(-p[:, 1])[:max_peaks]]
    p = p[np.argsort(p[:, 0])]
    return p.astype(np.float32)


# ----------------------------------------------------------------------------
# Molecular-formula enumeration from an accurate neutral mass (CHNOPS)
# ----------------------------------------------------------------------------
_EL = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956,
       "P": 30.97376163, "S": 31.972071}


def enumerate_formulas(mass, ppm=5.0, max_c=70, max_n=8, max_o=16, max_s=2, max_p=2, top=5):
    """Plausible CHNOPS formulas for a neutral monoisotopic mass (closed-shell, DBE >= 0).
    Returns a list of (formula, error_ppm) sorted by |error|."""
    tol = mass * ppm * 1e-6
    out = []
    for s in range(max_s + 1):
        for p in range(max_p + 1):
            for n in range(max_n + 1):
                for o in range(max_o + 1):
                    rest = mass - s * _EL["S"] - p * _EL["P"] - n * _EL["N"] - o * _EL["O"]
                    if rest < 12.0 - tol:
                        continue
                    for c in range(1, max_c + 1):
                        r = rest - c * 12.0
                        if r < -tol:
                            break
                        h = int(round(r / _EL["H"]))
                        if h < 0 or abs(r - h * _EL["H"]) > tol:
                            continue
                        dbe2 = 2 * c + 2 + n + p - h          # twice the double-bond equivalents
                        if dbe2 < 0 or dbe2 % 2 or h > 2 * c + n + p + 3:
                            continue
                        f = "".join(f"{el}{k if k > 1 else ''}" for el, k in
                                    (("C", c), ("H", h), ("N", n), ("O", o), ("P", p), ("S", s)) if k > 0)
                        out.append((f, (r - h * _EL["H"]) / mass * 1e6))
    out.sort(key=lambda t: abs(t[1]))
    return out[:top]


# ----------------------------------------------------------------------------
# Feature clustering and blank (background) matching for mzML data
# ----------------------------------------------------------------------------
def cluster_scans(scans, ppm=10.0):
    """Group MS2 scans with the same ion mode and precursor m/z (within ppm)."""
    scans = sorted(scans, key=lambda s: (s["mode"], s["precursor_mz"]))
    groups, cur = [], []
    for s in scans:
        if cur and (s["mode"] != cur[0]["mode"] or
                    s["precursor_mz"] - cur[0]["precursor_mz"] > cur[0]["precursor_mz"] * ppm * 1e-6):
            groups.append(cur)
            cur = []
        cur.append(s)
    if cur:
        groups.append(cur)
    return groups


def group_to_feature(group):
    """One representative scan (most peaks) + scan count and RT range of the group."""
    rep = dict(max(group, key=lambda s: len(s["peaks"])))
    rep["precursor_mz"] = float(np.median([s["precursor_mz"] for s in group]))
    rts = [s["rt"] for s in group if s.get("rt") is not None]
    rep["n_scans"] = len(group)
    rep["rt_range"] = f"{min(rts):.0f}-{max(rts):.0f}" if rts else "n/a"
    rep["rt_median"] = float(np.median(rts)) if rts else None
    rep["rt"] = rep.get("rt") if rep.get("rt") is not None else (rep["rt_median"] if rts else 0.0)
    return rep


def make_blank_matcher(blank_mz, ppm=10.0):
    arr = np.sort(np.asarray(list(blank_mz), dtype=np.float64))

    def in_blank(mz):
        if len(arr) == 0:
            return False
        tol = mz * ppm * 1e-6
        i = int(np.searchsorted(arr, mz))
        return any(0 <= j < len(arr) and abs(arr[j] - mz) <= tol for j in (i - 1, i))
    return in_blank


# ----------------------------------------------------------------------------
# mzML reading (pyopenms imported lazily so the other scripts do not need it)
# ----------------------------------------------------------------------------
class MzMLError(Exception):
    """A readable explanation of why a spectrum file could not be used."""


SUPPORTED_FORMATS = "mzML, mzXML, mzData, MGF, MSP (each optionally .gz compressed)"


def _check_mzml(path):
    """Cheap sanity checks before the file is parsed.
    Returns 'gzip' for compressed files, otherwise the format: mzml, mzxml, mzdata, mgf or msp."""
    import os
    import re
    size = os.path.getsize(path)
    if size == 0:
        raise MzMLError("the file is empty (0 bytes) - the copy or upload did not finish")
    with open(path, "rb") as f:
        head = f.read(8192)
        f.seek(max(size - 4096, 0))
        tail = f.read()
    if head[:2] == b"\x1f\x8b":
        return "gzip"
    if head[:3] == b"CDF":
        raise MzMLError("this is a netCDF (.cdf) file, which is not supported - convert it to mzML with ProteoWizard msconvert")
    if head[:4] == b"\x89HDF":
        raise MzMLError("this is an HDF5-based file (mzMLb or netCDF-4), which is not supported - convert it to "
                        "plain mzML with ProteoWizard msconvert")
    text = head.decode("utf-8", "ignore")
    end = tail.decode("utf-8", "ignore").rstrip()

    def closed(tag, *alt):
        if not any(end.endswith(t) for t in (tag,) + alt):
            raise MzMLError(f"the file is incomplete (it does not end with {tag}) - the transfer was probably cut "
                            "off; copy or upload it again")
    if "<mzML" in text or "<indexedmzML" in text:
        closed("</mzML>", "</indexedmzML>")
        return "mzml"
    if "<mzXML" in text:
        closed("</mzXML>")
        return "mzxml"
    if "<mzData" in text:
        closed("</mzData>")
        return "mzdata"
    if "BEGIN IONS" in text:
        return "mgf"
    if re.search(r"(?im)^\s*(name|num peaks)\s*:", text):
        return "msp"
    if text.lstrip().lower().startswith(("<!doctype html", "<html")):
        raise MzMLError("this is a web page (HTML), not a spectrum file - the download probably failed")
    raise MzMLError(f"this is not a recognised spectrum file. Supported: {SUPPORTED_FORMATS}. Raw vendor files "
                    "(.raw, .d, .wiff) must first be converted with ProteoWizard msconvert")


def _mzml_version(path):
    import re
    try:
        with open(path, "rb") as f:
            m = re.search(r'<mzML[^>]*?\sversion="([^"]+)"', f.read(8192).decode("utf-8", "ignore"))
        return m.group(1) if m else None
    except OSError:
        return None


def read_mzml(path, min_peaks=5, max_peaks=100, stats=None):
    """Read MS2 spectra from mzML, mzXML, mzData, MGF or MSP files (optionally .gz compressed).
    If a dict is passed as `stats` it is filled with scan counts (for helpful messages)."""
    import gzip
    import os
    import shutil
    import tempfile
    kind = _check_mzml(path)
    tmp = None
    n_ms2 = n_ms1 = n_noprec = n_few = 0
    out = []
    try:
        path_to_load = path
        if kind == "gzip":
            tmp = tempfile.NamedTemporaryFile(suffix=".dat", delete=False)
            with gzip.open(path, "rb") as src:
                shutil.copyfileobj(src, tmp)
            tmp.close()
            kind = _check_mzml(tmp.name)
            if kind == "gzip":
                raise MzMLError("the file is compressed twice")
            path_to_load = tmp.name

        if kind in ("mgf", "msp"):                                   # ---- peak lists
            with open(path_to_load, "r", encoding="utf-8", errors="replace") as fh:
                it = iter_mgf(fh) if kind == "mgf" else iter_msp(fh)
                for k, r in enumerate(it):
                    if r["ms_level"] != 2:
                        n_ms1 += r["ms_level"] == 1
                        continue
                    n_ms2 += 1
                    pmz = r["precursor_mz"]
                    if not pmz:
                        n_noprec += 1
                        continue
                    peaks = clean_peaks(r["peaks"], pmz, max_peaks)
                    if len(peaks) < min_peaks:
                        n_few += 1
                        continue
                    out.append({"scan_index": k, "rt": r.get("rt"), "precursor_mz": float(pmz), "peaks": peaks,
                                "mode": r.get("mode") or "Positive", "ce": r.get("ce")})
        else:                                                        # ---- XML formats via OpenMS
            import pyopenms as oms
            exp = oms.MSExperiment()
            loader = {"mzxml": oms.MzXMLFile, "mzdata": oms.MzDataFile}.get(kind, oms.MzMLFile)()
            try:
                loader.load(path_to_load, exp)
            except Exception as exc:
                first = str(exc).strip().splitlines()[0][:150] if str(exc).strip() else type(exc).__name__
                ver = _mzml_version(path_to_load) if kind == "mzml" else None
                if ver and ver.startswith("0."):
                    raise MzMLError(f"the file declares mzML version {ver}, an old pre-release format that OpenMS "
                                    "cannot read (typical for test files shipped with software packages). Use a "
                                    "real data file, or convert it with ProteoWizard msconvert") from None
                raise MzMLError(f"OpenMS could not parse the file ({first}). It may be corrupted or an unusual "
                                f"{kind} variant - try re-converting it with msconvert (mzML, 32-bit, zlib "
                                "compression)") from None
            for i, spec in enumerate(exp):
                if spec.getMSLevel() != 2:
                    n_ms1 += spec.getMSLevel() == 1
                    continue
                n_ms2 += 1
                precs = spec.getPrecursors()
                if not precs:
                    n_noprec += 1
                    continue
                pmz = precs[0].getMZ()
                mz, inten = spec.get_peaks()
                peaks = clean_peaks(np.stack([mz, inten], 1), pmz, max_peaks)
                if len(peaks) < min_peaks:
                    n_few += 1
                    continue
                pol = str(spec.getInstrumentSettings().getPolarity()).upper()
                mode = "Negative" if ("NEG" in pol or pol.strip() == "2") else "Positive"
                ce = precs[0].getActivationEnergy() or None
                out.append({"scan_index": i, "rt": spec.getRT(), "precursor_mz": pmz,
                            "peaks": peaks, "mode": mode, "ce": ce})
    finally:
        if tmp is not None:
            try:
                os.unlink(tmp.name)
            except OSError:
                pass
    print(f"{path}: {n_ms2} MS2 spectra found ({kind})")
    if stats is not None:
        stats.update(ms1=n_ms1, ms2=n_ms2, no_precursor=n_noprec, too_few_peaks=n_few, kept=len(out), format=kind)
    return out


# ----------------------------------------------------------------------------
# Spectral library search (cosine) against the cleaned dataset
# ----------------------------------------------------------------------------
def cosine_score(q, l, tol=0.02):
    """Greedy one-to-one cosine between two (n,2) peak arrays (sqrt intensities).
    Returns (score, n_matched)."""
    qm, li = q[:, 0], l[:, 0]
    d = np.abs(qm[:, None] - li[None, :])
    qi, lj = np.where(d <= tol)
    if len(qi) == 0:
        return 0.0, 0
    qs, ls = np.sqrt(q[:, 1]), np.sqrt(l[:, 1])
    prod = qs[qi] * ls[lj]
    order = np.argsort(-prod)
    used_q, used_l, tot, n = set(), set(), 0.0, 0
    for k in order:
        a, b = qi[k], lj[k]
        if a in used_q or b in used_l:
            continue
        used_q.add(a); used_l.add(b)
        tot += prod[k]; n += 1
    denom = np.sqrt((qs ** 2).sum() * (ls ** 2).sum())
    return float(tot / denom) if denom > 0 else 0.0, n


class LibraryIndex:
    def __init__(self, records):
        self.recs = sorted(records, key=lambda r: (r["mode"], r["precursor_mz"]))
        self.mode = np.array([r["mode"] for r in self.recs])
        self.mz = np.array([r["precursor_mz"] for r in self.recs])

    def query(self, feat, ppm=20.0, frag_tol=0.02, min_matched=3, top=5, min_score=0.4):
        sel = np.where(self.mode == feat["mode"])[0]
        if len(sel) == 0:
            return [], 0
        mzs = self.mz[sel]
        tol = feat["precursor_mz"] * ppm * 1e-6
        lo, hi = np.searchsorted(mzs, [feat["precursor_mz"] - tol, feat["precursor_mz"] + tol])
        hits = []
        for k in sel[lo:hi]:
            r = self.recs[k]
            s, n = cosine_score(feat["peaks"], r["peaks"], frag_tol)
            if n >= min_matched and s >= min_score:
                hits.append((s, n, r))
        hits.sort(key=lambda t: -t[0])
        # one hit per molecule
        seen, out = set(), []
        for s, n, r in hits:
            if r["inchikey14"] in seen:
                continue
            seen.add(r["inchikey14"])
            out.append((s, n, r))
            if len(out) >= top:
                break
        return out, hi - lo


# ----------------------------------------------------------------------------
# Retrieval pool: candidate structures searched by mass, ranked by fingerprint
# ----------------------------------------------------------------------------
class RetrievalPool:
    """items: iterable of (smiles, inchikey14, exact_mass). One entry per inchikey14."""

    def __init__(self, items, verbose=False):
        seen, rows = set(), []
        for smi, key, mass in items:
            if key in seen:
                continue
            seen.add(key)
            rows.append((mass, key, smi))
        rows.sort()
        self.mass = np.array([r[0] for r in rows])
        self.key = [r[1] for r in rows]
        self.smiles = [r[2] for r in rows]
        bits = []
        for i, (_, _, smi) in enumerate(rows):
            bits.append(np.packbits(morgan_bits(Chem.MolFromSmiles(smi))))
            if verbose and (i + 1) % 50000 == 0:
                print(f"    fingerprints {i + 1:,}/{len(rows):,}", flush=True)
        self.bits = np.stack(bits) if bits else np.zeros((0, 256), np.uint8)

    def __len__(self):
        return len(self.key)

    @classmethod
    def from_records(cls, recs, verbose=False):
        return cls(((r["smiles"], r["inchikey14"], r["exact_mass"]) for r in recs), verbose)

    def rank(self, p, neutral_mass_, ppm, true_key):
        """p: predicted fingerprint probabilities (2048,). Returns (n_candidates, rank of truth or None)."""
        tol = neutral_mass_ * ppm * 1e-6
        lo, hi = np.searchsorted(self.mass, [neutral_mass_ - tol, neutral_mass_ + tol])
        if hi <= lo:
            return 0, None
        cand = np.unpackbits(self.bits[lo:hi], axis=1).astype(np.float32)
        inter = cand @ p
        score = inter / (cand.sum(1) + p.sum() - inter + 1e-9)
        order = np.argsort(-score)
        keys = [self.key[lo + j] for j in order]
        try:
            return hi - lo, keys.index(true_key) + 1
        except ValueError:
            return hi - lo, None


def _pool_top_candidates(self, p, neutral_mass_, ppm, top=10):
    """Best `top` candidates inside the mass window, ranked by soft Tanimoto to the predicted fingerprint."""
    tol = neutral_mass_ * ppm * 1e-6
    lo, hi = np.searchsorted(self.mass, [neutral_mass_ - tol, neutral_mass_ + tol])
    if hi <= lo:
        return [], 0
    cand = np.unpackbits(self.bits[lo:hi], axis=1).astype(np.float32)
    inter = cand @ p
    score = inter / (cand.sum(1) + p.sum() - inter + 1e-9)
    order = np.argsort(-score)[:top]
    return [{"smiles": self.smiles[lo + j], "inchikey14": self.key[lo + j],
             "exact_mass": float(self.mass[lo + j]), "score": float(score[j])} for j in order], int(hi - lo)


def _pool_save(self, path):
    import pickle
    with open(path, "wb") as f:
        pickle.dump({"mass": self.mass, "key": self.key, "smiles": self.smiles, "bits": self.bits}, f, protocol=4)


def _pool_load(cls, path):
    import pickle
    with open(path, "rb") as f:
        d = pickle.load(f)
    obj = cls.__new__(cls)
    obj.mass, obj.key, obj.smiles, obj.bits = d["mass"], d["key"], d["smiles"], d["bits"]
    return obj


RetrievalPool.top_candidates = _pool_top_candidates
RetrievalPool.save = _pool_save
RetrievalPool.load = classmethod(_pool_load)


def retrieval_metrics(results, ks=(1, 5, 10, 50)):
    """results: list of (n_candidates, rank or None) -> dict incl. random baseline and per-size breakdown."""
    nq = len(results)
    out = {"queries": nq}
    for k in ks:
        out[f"top{k}"] = sum(1 for n, r in results if r is not None and r <= k) / max(nq, 1)
        out[f"random_top{k}"] = sum(min(k, n) / n for n, r in results if n > 0) / max(nq, 1)
    ns = [n for n, _ in results]
    out["median_candidates"] = float(np.median(ns)) if ns else 0
    out["mean_candidates"] = float(np.mean(ns)) if ns else 0
    out["truth_outside_mass_window"] = sum(1 for n, r in results if r is None) / max(nq, 1)
    bins = [(1, 1), (2, 10), (11, 100), (101, 1000), (1001, 10**9)]
    by = {}
    for lo, hi in bins:
        sub = [(n, r) for n, r in results if lo <= n <= hi]
        if sub:
            by[f"{lo}-{hi if hi < 10**9 else 'inf'}"] = {
                "queries": len(sub),
                "top1": sum(1 for n, r in sub if r == 1) / len(sub),
                "top10": sum(1 for n, r in sub if r is not None and r <= 10) / len(sub),
                "random_top1": sum(1 / n for n, r in sub) / len(sub)}
    out["by_candidate_count"] = by
    return out


# ----------------------------------------------------------------------------
# Adduct inference from co-eluting precursors
# ----------------------------------------------------------------------------
ADDUCTS_BY_MODE = {
    "Positive": ["[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+"],
    "Negative": ["[M-H]-", "[M+FA-H]-", "[M+Cl]-"],
}


def infer_adducts(feats, ppm=5.0, rt_tol=10.0):
    """feats: dicts with precursor_mz, mode, rt_median. Two co-eluting precursors whose neutral masses agree
    under two different adducts (e.g. [M+H]+ / [M+NH4]+) are assigned those adducts.
    Returns a list of (adduct, evidence) per feature; unpaired features get the default adduct, flagged 'assumed'."""
    n = len(feats)
    out = [None] * n
    for i in range(n):
        for j in range(i + 1, n):
            fi, fj = feats[i], feats[j]
            if fi["rt_median"] is None or fj["rt_median"] is None:
                continue                      # without retention times a pairing would be guesswork
            if fi["mode"] != fj["mode"] or abs(fi["rt_median"] - fj["rt_median"]) > rt_tol:
                continue
            best = None
            for a1 in ADDUCTS_BY_MODE.get(fi["mode"], []):
                for a2 in ADDUCTS_BY_MODE.get(fi["mode"], []):
                    if a1 == a2:
                        continue
                    m1 = fi["precursor_mz"] - ADDUCT_DELTA[a1]
                    m2 = fj["precursor_mz"] - ADDUCT_DELTA[a2]
                    err = abs(m1 - m2) / m1 * 1e6
                    if err <= ppm and (best is None or err < best[0]):
                        best = (err, a1, a2)
            if best:
                err, a1, a2 = best
                if out[i] is None:
                    out[i] = (a1, f"co-eluting pair with m/z {fj['precursor_mz']:.4f} ({a2}), {err:.1f} ppm")
                if out[j] is None:
                    out[j] = (a2, f"co-eluting pair with m/z {fi['precursor_mz']:.4f} ({a1}), {err:.1f} ppm")
    for i in range(n):
        if out[i] is None:
            default = "[M-H]-" if feats[i]["mode"] == "Negative" else "[M+H]+"
            out[i] = (default, "assumed")
    return out


def confidence_label(score, n_matched):
    if score >= 0.7 and n_matched >= 6:
        return "strong"
    if score >= 0.5 and n_matched >= 4:
        return "moderate"
    return "weak"
