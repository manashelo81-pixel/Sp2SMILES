"""Look inside a MoNA export (.msp, or a .zip containing one) WITHOUT extracting it.
Prints the raw first record, which fields exist, and how many records have usable
structures / adducts / ion modes, so you can check the parser fits your export.

  python 00_inspect_mona.py \
     --input /home/manas/Desktop/SMILES/experimental/MoNA-export-LC-MS_Spectra-msp.zip --n 20000
"""
import argparse
import io
import zipfile
from collections import Counter

from rdkit import Chem, RDLogger

import common as C

RDLogger.DisableLog("rdApp.*")


def open_text(path):
    if path.lower().endswith(".zip"):
        z = zipfile.ZipFile(path)
        name = next(n for n in z.namelist() if n.lower().endswith(".msp"))
        print("zip member:", name)
        return io.TextIOWrapper(z.open(name), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--n", type=int, default=20000, help="records to scan")
    args = ap.parse_args()

    # ---- raw view: first record + field-name frequencies ----
    keys, first, nrec, cur = Counter(), [], 0, []
    with open_text(args.input) as f:
        for line in f:
            s = line.rstrip("\n")
            if not s.strip():
                if cur:
                    nrec += 1
                    if nrec == 1:
                        first = cur[:40]
                    for l in cur:
                        if ":" in l and not l[0].isdigit():
                            keys[l.split(":", 1)[0].strip()] += 1
                    cur = []
                    if nrec >= args.n:
                        break
            else:
                cur.append(s)
    print("\n===== FIRST RECORD (first 40 lines) =====")
    print("\n".join(l[:200] for l in first))
    print(f"\n===== FIELDS seen in {nrec:,} records =====")
    for k, v in keys.most_common(40):
        print(f"{k:28s} {v:8,d}  ({v / max(nrec, 1):.0%})")

    # ---- parsed view (what 01_prepare_data.py will actually use) ----
    st = Counter()
    modes, adducts, inst, cls = Counter(), Counter(), Counter(), Counter()
    n = 0
    with open_text(args.input) as f:
        for r in C.iter_msp(f):
            n += 1
            st["has_smiles"] += bool(r["smiles_raw"])
            st["has_inchi"] += bool(r["inchi_raw"])
            st["has_ce"] += r["ce"] is not None
            st["ms2"] += r["ms_level"] == 2
            ok = bool(r["smiles_raw"]) and Chem.MolFromSmiles(r["smiles_raw"]) is not None
            ok2 = ok or (r["inchi_raw"].startswith("InChI=") and Chem.MolFromInchi(r["inchi_raw"]) is not None)
            st["structure_usable"] += ok2
            st["precursor_mz"] += bool(r["precursor_mz"])
            modes[r["mode"]] += 1
            adducts[r["adduct_raw"] or "(none)"] += 1
            inst[r["instrument"]] += 1
            if n >= min(args.n, 20000):
                break
    print(f"\n===== PARSED (first {n:,}) =====")
    for k, v in st.items():
        print(f"{k:18s} {v:8,d}  ({v / max(n, 1):.0%})")
    print("ion modes:", dict(modes))
    print("precursor types:", dict(adducts.most_common(10)))
    print("instrument classes:", dict(inst))
    print("\nIf 'structure_usable' is low, the SMILES/InChI are stored somewhere the parser does not look;"
          " paste the FIRST RECORD output to me and I will adapt iter_msp().")


if __name__ == "__main__":
    main()
