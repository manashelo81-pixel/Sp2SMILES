"""Predict SMILES for MS/MS spectra in mzML files using a trained model.

  python 03_predict.py --ckpt /home/manas/Desktop/SMILES/run1/model.pt \
      --mzml /home/manas/Desktop/SMILES/A1_GA1_01_2600.mzML \
      --out predictions.csv
Use --blank blank1.mzML [blank2.mzML ...] to flag background ions seen in your blanks.
Optional: --candidates my_candidates.smi  (one SMILES per line, e.g. exported from PubChem)
adds database-style fingerprint retrieval on top of generation.
"""
import argparse
import csv

import numpy as np
import torch
from rdkit import Chem

import common as C
import model as M


read_mzml = C.read_mzml


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--mzml", nargs="+", required=True)
    ap.add_argument("--out", default="predictions.csv")
    ap.add_argument("--n-samples", type=int, default=64)
    ap.add_argument("--temperature", type=float, default=1.0)
    ap.add_argument("--ppm", type=float, default=10.0, help="mass tolerance for candidates (ppm)")
    ap.add_argument("--mass-tol", type=float, default=None, help="tolerance in Da (overrides --ppm)")
    ap.add_argument("--formula-ppm", type=float, default=5.0, help="tolerance for the formula list")
    ap.add_argument("--allow-charged", action="store_true", help="keep net-charged generated structures")
    ap.add_argument("--top", type=int, default=10, help="candidates reported per spectrum")
    ap.add_argument("--min-peaks", type=int, default=5)
    ap.add_argument("--candidates", help="optional file with one SMILES per line")
    ap.add_argument("--max-spectra", type=int, default=0, help="limit number of unique precursors (testing)")
    ap.add_argument("--blank", nargs="+", default=[],
                    help="blank mzML file(s): precursors also present in a blank are flagged as background and skipped")
    ap.add_argument("--blank-ppm", type=float, default=10.0)
    ap.add_argument("--cluster-ppm", type=float, default=10.0,
                    help="scans with the same precursor (within ppm) are predicted once")
    ap.add_argument("--no-cluster", action="store_true", help="predict every scan separately")
    ap.add_argument("--instrument", default="UNKNOWN", help="ORBITRAP, QTOF, QQQ, ITRAP, OTHER or UNKNOWN")
    ap.add_argument("--collision-energy", type=float, default=None,
                    help="overrides the value read from the mzML (if any)")
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    net, ck = M.load_model(args.ckpt, device)
    itos = {i: t for i, t in enumerate(ck["vocab"])}
    ad = {a: i for i, a in enumerate(ck["adducts"])}
    ins = {a: i for i, a in enumerate(ck.get("instruments", ["UNKNOWN"]))}

    pool = []
    if args.candidates:
        seen = set()
        for line in open(args.candidates):
            mol, flat = C.standardize(line.split()[0] if line.strip() else "")
            if mol is not None and flat not in seen:
                seen.add(flat)
                info = C.mol_info(mol)
                pool.append((flat, mol, info["exact_mass"], info["inchikey14"], info["formula"],
                             C.morgan_bits(mol).astype(np.float32)))
        pool.sort(key=lambda x: x[2])
        pool_mass = np.array([p[2] for p in pool])
        print(f"{len(pool)} candidate structures loaded")

    blank_mz = []
    for bpath in args.blank:
        blank_mz += [s["precursor_mz"] for s in read_mzml(bpath, 1, 100)]
    in_blank = C.make_blank_matcher(blank_mz, args.blank_ppm)
    if args.blank:
        print(f"{len(blank_mz)} blank MS2 scans loaded from {len(args.blank)} file(s)")

    rows = []
    fcache = {}
    done = 0
    n_bg = 0
    for path in args.mzml:
        scans = list(read_mzml(path, args.min_peaks, 100))
        groups = [[s] for s in scans] if args.no_cluster else C.cluster_scans(scans, args.cluster_ppm)
        feats = [C.group_to_feature(g) for g in groups]
        print(f"{path}: {len(scans)} scans -> {len(feats)} unique precursors")
        for s in feats:
            if args.max_spectra and done >= args.max_spectra:
                break
            done += 1
            adduct = "[M-H]-" if s["mode"] == "Negative" else "[M+H]+"
            if in_blank(s["precursor_mz"]):
                n_bg += 1
                rows.append([path, s["scan_index"], round(s["rt"], 2), round(s["precursor_mz"], 5), adduct,
                             0, "", "", "", "", 0, "", "background: precursor also found in blank", "",
                             s["n_scans"], s["rt_range"]])
                continue
            nm = C.neutral_mass(s["precursor_mz"], adduct)
            tol = args.mass_tol if args.mass_tol is not None else nm * args.ppm * 1e-6
            fkey = round(nm, 4)
            if fkey not in fcache:
                fcache[fkey] = ";".join(f"{f}({e:+.1f}ppm)" for f, e in
                                        C.enumerate_formulas(nm, args.formula_ppm))
            fcands = fcache[fkey]
            ce = args.collision_energy if args.collision_energy is not None else s["ce"]
            x = M.make_batch([(s["peaks"], s["precursor_mz"], ad.get(adduct, 0),
                               ins.get(args.instrument.upper(), 0), ce)], device)
            with torch.no_grad():
                mem, kp, fpl = net.encode(*x)
                p = torch.sigmoid(fpl)[0].cpu().numpy()
                greedy = M.generate(net, mem, kp, ck["max_tokens"] + 5, 0.0)
                mem_n, kp_n = mem.repeat(args.n_samples, 1, 1), kp.repeat(args.n_samples, 1)
                samp = M.generate(net, mem_n, kp_n, ck["max_tokens"] + 5, args.temperature)
            cands = {}
            for row in list(greedy) + list(samp):
                flat, mol = M.tokens_to_mol(row, itos)
                if mol is None:
                    continue
                if not args.allow_charged and Chem.GetFormalCharge(mol) != 0:
                    continue
                info = C.mol_info(mol)
                err = info["exact_mass"] - nm
                if abs(err) > tol:
                    continue
                c = cands.setdefault(flat, {"mol": mol, "info": info, "n": 0, "src": "generated"})
                c["n"] += 1
            if pool:
                lo, hi = np.searchsorted(pool_mass, [nm - tol, nm + tol])
                for flat, mol, em, ik, fm, bits in pool[lo:hi]:
                    cands.setdefault(flat, {"mol": mol, "info": {"exact_mass": em, "formula": fm,
                                                                  "inchikey14": ik}, "n": 0, "src": "database"})
            scored = []
            for flat, c in cands.items():
                bits = C.morgan_bits(c["mol"]).astype(np.float32)
                inter = float(bits @ p)
                fps = inter / (bits.sum() + p.sum() - inter + 1e-9)
                scored.append((c["n"], fps, flat, c))
            scored.sort(key=lambda t: (-t[0], -t[1]))
            if not scored:
                rows.append([path, s["scan_index"], round(s["rt"], 2), round(s["precursor_mz"], 5),
                             adduct, 0, "", "", "", "", 0, "", "no mass-consistent candidate", fcands,
                             s["n_scans"], s["rt_range"]])
            for rank, (n, fps, flat, c) in enumerate(scored[:args.top], 1):
                ppm = (c["info"]["exact_mass"] - nm) / nm * 1e6
                rows.append([path, s["scan_index"], round(s["rt"], 2), round(s["precursor_mz"], 5),
                             adduct, rank, flat, c["info"]["formula"], round(ppm, 2),
                             c["src"], n, round(fps, 3), "", fcands, s["n_scans"], s["rt_range"]])
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "scan_index", "rt_s", "precursor_mz", "adduct", "rank", "smiles",
                    "formula", "mass_error_ppm", "source", "times_generated", "fp_score", "note", "possible_formulas", "n_scans", "rt_range_s"])
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows for {done} unique precursors ({n_bg} flagged as background) -> {args.out}")


if __name__ == "__main__":
    main()
