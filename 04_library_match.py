"""Library (spectral) search of your mzML features against the cleaned MoNA/GNPS dataset.
Known compounds that are in the library can be identified directly, with no model needed.

  python 04_library_match.py --dataset /home/manas/Desktop/SMILES/processed_all/dataset.pkl \
      --mzml /home/manas/Desktop/SMILES/A1_GA1_01_2600.mzML \
      --blank /home/manas/Desktop/SMILES/BE8-blank_BE8_01_2612.mzML /home/manas/Desktop/SMILES/BE8-blank_BE8_01_2625.mzML \
      --out library_matches.csv
"""
import argparse
import csv
import pickle

import common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, help="dataset.pkl from 01_prepare_data.py")
    ap.add_argument("--mzml", nargs="+", required=True)
    ap.add_argument("--blank", nargs="+", default=[])
    ap.add_argument("--out", default="library_matches.csv")
    ap.add_argument("--ppm", type=float, default=20.0, help="precursor tolerance (ppm)")
    ap.add_argument("--frag-tol", type=float, default=0.02, help="fragment tolerance (Da)")
    ap.add_argument("--min-matched", type=int, default=3)
    ap.add_argument("--min-score", type=float, default=0.4)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--blank-ppm", type=float, default=10.0)
    ap.add_argument("--cluster-ppm", type=float, default=10.0)
    args = ap.parse_args()

    print("loading library ...", flush=True)
    recs = pickle.load(open(args.dataset, "rb"))
    lib = C.LibraryIndex(recs)
    print(f"{len(recs):,} library spectra")

    blank_mz = []
    for b in args.blank:
        blank_mz += [s["precursor_mz"] for s in C.read_mzml(b, 1, 100)]
    in_blank = C.make_blank_matcher(blank_mz, args.blank_ppm)

    rows, n_feat, n_bg, n_hit = [], 0, 0, 0
    for path in args.mzml:
        scans = C.read_mzml(path, 3, 100)
        feats = [C.group_to_feature(g) for g in C.cluster_scans(scans, args.cluster_ppm)]
        print(f"{path}: {len(scans)} scans -> {len(feats)} unique precursors")
        for f in feats:
            n_feat += 1
            base = [path.split("/")[-1], round(f["precursor_mz"], 5), f["rt_range"], f["n_scans"]]
            if in_blank(f["precursor_mz"]):
                n_bg += 1
                rows.append(base + [0, "", "", "", "", "", "", "", "", "background: also in blank"])
                continue
            hits, n_cand = lib.query(f, args.ppm, args.frag_tol, args.min_matched, args.top, args.min_score)
            if not hits:
                rows.append(base + [0, "", "", "", "", "", "", "", "",
                                    f"no match ({n_cand} library spectra at this precursor mass)"])
                continue
            n_hit += 1
            for rank, (s, n, r) in enumerate(hits, 1):
                ppm = (r["precursor_mz"] - f["precursor_mz"]) / f["precursor_mz"] * 1e6
                rows.append(base + [rank, round(s, 3), n, r["smiles"], r["formula"], r["adduct"],
                                    r["instrument"], r["id"], round(ppm, 1), ""])
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "precursor_mz", "rt_range_s", "n_scans", "rank", "cosine", "matched_peaks",
                    "smiles", "formula", "library_adduct", "library_instrument", "library_id",
                    "precursor_diff_ppm", "note"])
        w.writerows(rows)
    print(f"{n_feat} unique precursors: {n_bg} background, {n_hit} with library matches -> {args.out}")


if __name__ == "__main__":
    main()
