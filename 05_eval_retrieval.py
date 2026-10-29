"""Realistic evaluation of the fingerprint-retrieval model on held-out (unseen-scaffold) spectra.

For each test spectrum the model predicts a fingerprint; candidate structures are every molecule in the pool whose
mass fits the precursor (within --ppm); candidates are ranked by soft Tanimoto to the prediction.
Pools: 'test' (only test molecules - optimistic), 'all' (every library molecule - more realistic) and optionally
'file' (your own candidate list, e.g. PubChem structures; the true molecule is added if missing).

  python 05_eval_retrieval.py --ckpt /home/manas/Desktop/SMILES/run2/model.pt \
      --dataset /home/manas/Desktop/SMILES/processed_all/dataset.pkl --out retrieval_eval.json
"""
import argparse
import json
import pickle

import numpy as np
import torch

import common as C
import model as M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--split", default="test")
    ap.add_argument("--max-queries", type=int, default=5000)
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--pool-file", help="optional file with one SMILES per line")
    ap.add_argument("--out", default="retrieval_eval.json")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    net, ck = M.load_model(args.ckpt, device)
    ad = {a: i for i, a in enumerate(ck["adducts"])}
    ins = {a: i for i, a in enumerate(ck.get("instruments", ["UNKNOWN"]))}

    recs = pickle.load(open(args.dataset, "rb"))
    queries = [r for r in recs if r["split"] == args.split]
    rng = np.random.default_rng(args.seed)
    if len(queries) > args.max_queries:
        queries = [queries[i] for i in rng.choice(len(queries), args.max_queries, replace=False)]
    print(f"{len(queries):,} {args.split} queries from {len(recs):,} library spectra")

    # ---- model: predicted fingerprints ----
    probs = []
    net.eval()
    with torch.no_grad():
        for s in range(0, len(queries), 128):
            chunk = queries[s:s + 128]
            x = M.make_batch([(r["peaks"], r["precursor_mz"], ad.get(r["adduct"], 0),
                               ins.get(r.get("instrument", "UNKNOWN"), 0), r.get("ce")) for r in chunk], device)
            _, _, fpl = net.encode(*x)
            probs.append(torch.sigmoid(fpl).float().cpu().numpy())
    probs = np.concatenate(probs)

    # ---- pools ----
    pools = {}
    print("building pool: test molecules ...", flush=True)
    pools["test_molecules_only"] = C.RetrievalPool.from_records([r for r in recs if r["split"] == args.split])
    print("building pool: all library molecules ...", flush=True)
    pools["all_library_molecules"] = C.RetrievalPool.from_records(recs, verbose=True)
    if args.pool_file:
        from rdkit import Chem
        items = []
        for line in open(args.pool_file):
            mol, flat = C.standardize(line.split()[0] if line.strip() else "")
            if mol is not None:
                info = C.mol_info(mol)
                items.append((flat, info["inchikey14"], info["exact_mass"]))
        items += [(r["smiles"], r["inchikey14"], r["exact_mass"]) for r in queries]   # truth always present
        print(f"building pool: file ({len(items):,} structures incl. truths) ...", flush=True)
        pools["pool_file_plus_truths"] = C.RetrievalPool(items, verbose=True)

    report = {"split": args.split, "ppm": args.ppm, "n_library_spectra": len(recs)}
    for name, pool in pools.items():
        res = []
        for r, p in zip(queries, probs):
            nm = C.neutral_mass(r["precursor_mz"], r["adduct"])
            res.append(pool.rank(p, nm, args.ppm, r["inchikey14"]))
        m = C.retrieval_metrics(res)
        m["pool_size"] = len(pool)
        report[name] = m
        print(f"\n=== pool: {name} ({len(pool):,} molecules) ===")
        print(f"candidates per spectrum: median {m['median_candidates']:.0f}, mean {m['mean_candidates']:.1f}")
        for k in (1, 5, 10, 50):
            print(f"  top-{k:<3d} model {m[f'top{k}']:.3f}   random guessing {m[f'random_top{k}']:.3f}")
        print("  by number of candidates:")
        for b, v in m["by_candidate_count"].items():
            print(f"    {b:>9s}: n={v['queries']:5d}  top1 {v['top1']:.3f} (random {v['random_top1']:.3f})  top10 {v['top10']:.3f}")
    json.dump(report, open(args.out, "w"), indent=2)
    print(f"\nsaved {args.out}")


if __name__ == "__main__":
    main()
